package main 

import "core:fmt"
import "core:math"
import stbi "vendor:stb/image"

make_2D :: proc($T: typeid, M, N: int, allocator := context.allocator) -> [][]T {
	backing := make([]T, M*N)
	grid := make([][]T, N)
	for j  in 0..<N {
		grid[j] = backing[j*M:(j+1)*M]
	}
	return grid
}

delete_2D :: proc(grid: [][]$T) {
	delete(grid[0])
	delete(grid)
}

write_velocity :: proc(velocity: [][][2]f32, filename: cstring, SS: int) {
	N := len(velocity)
	assert(len(velocity[0]) == N)
	velocity_image := make_2D([4]u8, SS*N, SS*N)
	defer delete_2D(velocity_image)

	for j in 0..<N {
		for i in 0..<N {
			pixel_value := [4]u8{
				// Flowmap color map
				u8(255*(0.5 + 0.5*velocity[j][i].x)),
				u8(255*(0.5 + 0.5*velocity[j][i].y)),
				0,
				255,
			}

			// Supersample to make the image more legible
			for y in 0..<SS {
				for x in 0..<SS {
					velocity_image[SS*j+y][SS*i+x] = pixel_value
				}
			}
		}
	}
	stbi.write_png(filename, i32(SS*N), i32(SS*N), 4, &velocity_image[0][0], 4*i32(SS*N))
}

logarithmic_scalar_color_map :: proc(v: f32) -> [4]u8 {
    mix :: proc(a, b, x: [3]f32) -> [3]f32 {
    	return a * (1.0 - x) + b * x
    }
    logv := math.log(abs(v), 10.0)

    f := math.floor(logv + 7.0);
    i := math.floor(16 * ((logv + 7.0) - f));

    c: [4]f32 = {0, 0, 0, 1}
         if f < 0.0 do c.rgb = {0.0, 0.0, 0.0};
    else if f < 1.0 do c.rgb = mix({1.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, i / 16.0);
    else if f < 2.0 do c.rgb = mix({0.0, 1.0, 0.0}, {1.0, 1.0, 1.0}, i / 16.0);
    else if f < 3.0 do c.rgb = mix({0.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, i / 16.0);
    else if f < 4.0 do c.rgb = mix({1.0, 1.0, 0.0}, {1.0, 1.0, 1.0}, i / 16.0);
    else if f < 5.0 do c.rgb = mix({1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, i / 16.0);
    else if f < 6.0 do c.rgb = mix({0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, i / 16.0);
    else if f < 7.0 do c.rgb = mix({1.0, 0.5, 0.0}, {1.0, 1.0, 1.0}, i / 16.0);
    else if f < 8.0 do c.rgb = mix({1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, i / 16.0);

    return {
    	u8(255*c.r),
    	u8(255*c.g),
    	u8(255*c.b),
    	u8(255*c.a),
    }
}

linear_scalar_color_map :: proc(f: f32) -> [4]u8 {
	return [4]u8{
		u8(clamp(-255*f, 0, 255)), // negative values are red
		u8(clamp(+255*f, 0, 255)), // positive values are green
		0,
		255,
	}
}

write_scalar :: proc(scalar: [][]f32, filename: cstring, SS: int, color_map: proc(f32) -> [4]u8, scale: f32 = 1.0) {
	N := len(scalar)
	assert(len(scalar[0]) == N)
	scalar_image := make_2D([4]u8, SS*N, SS*N)
	defer delete_2D(scalar_image)

	max_value: f32 = 0.0
	if scale != 1.0 {
		for j in 0..<N {
			for i in 0..<N {
				max_value = max(max_value, abs(scalar[j][i]))
			}
		}
	} else {
		max_value = 1.0
	}

	for j in 0..<N {
		for i in 0..<N {
			pixel_value := color_map((1.0/max_value)*scalar[j][i])

			// Supersample to make the image more legible
			for y in 0..<SS {
				for x in 0..<SS {
					scalar_image[SS*j+y][SS*i+x] = pixel_value
				}
			}
		}
	}
	stbi.write_png(filename, i32(SS*N), i32(SS*N), 4, &scalar_image[0][0], 4*i32(SS*N))
}

MAIN_N :: 32

main :: proc() {
	stbi.flip_vertically_on_write(true)

	/*
	finest := make_2D(f32, 8, 8)
	fine := make_2D(f32, 4, 4)
	coarse := make_2D(f32, 2, 2)
	coarse[1][1] = 1.0
	prolongate(coarse, fine)
	for j in 0..<4 do fmt.println(fine[j])
	prolongate(fine, finest)
	for j in 0..<8 do fmt.println(finest[j])
	if true do return
	*/
	N := MAIN_N
	dump_interval := N
	iterations := N*N
	iterations = 1
	assert(N & (N-1) == 0, "N not power of two")

	// Initialize Velocity
	initial_velocity := make_2D([2]f32, N, N)
	defer delete_2D(initial_velocity)
	for j in 3*N/8..<5*N/8 {
		for i in 3*N/8..<5*N/8 {
			initial_velocity[j][i].y = +1.0
		}
	}
	write_velocity(initial_velocity, "initial_velocity.png", 384/N)

	// The pressure and divergence grids are (N+1)x(N+1) in the vertex grid with 
	// NxN velocity cells. We still store these in NxN arrays, where we explicitly
	// maintain the lower boundary value and implicitly assume the upper boundary.
	initial_divergence := make_2D(f32, N, N)
	pressure_ping := make_2D(f32, N, N)
	pressure_pong := make_2D(f32, N, N)
	defer {
		delete_2D(initial_divergence)
		delete_2D(pressure_ping)
		delete_2D(pressure_pong)
	}

	calc_divergence(initial_velocity, initial_divergence)
	write_scalar(initial_divergence, "initial_divergence.png", 384/N, linear_scalar_color_map)

	residual0 := make_2D(f32, N, N)
	defer delete_2D(residual0)	
	calc_residual(pressure_ping, initial_divergence, residual0)	
	max_residual0, rms_residual0, gradient_residual0 := calc_stats(residual0)


	// Solve pressure
	for iter := 0; iter < ITER_MAX; iter += 1 {
	//for iter := 0; iter < 4; iter += 1 {
		fmt.println(iter)
		//iterate(pressure_ping, initial_divergence, pressure_pong, 4.0/5.0)
		//pressure_ping, pressure_pong = pressure_pong, pressure_ping
		vcycle(&pressure_ping, &initial_divergence, &pressure_pong, 0, N*N, 0)

		if !((iter % dump_interval) == 0) do continue
		write_scalar(pressure_ping, "pressure.png", 384/N, logarithmic_scalar_color_map)
		write_scalar(pressure_ping, "pressure_linear.png", 384/N, linear_scalar_color_map, 2.0)

		residual := make_2D(f32, N, N)
		defer delete_2D(residual)	
		calc_residual(pressure_ping, initial_divergence, residual)	
		write_scalar(residual, "residual.png", 384/N, logarithmic_scalar_color_map)

		velocity := make_2D([2]f32, N, N)
		defer delete_2D(velocity)
		subtract_gradient(pressure_ping, initial_velocity, velocity)
		write_velocity(velocity, "velocity.png", 384/N)

		divergence := make_2D(f32, N, N)
		defer delete_2D(divergence)
		calc_divergence(velocity, divergence)
		write_scalar(divergence, "divergence.png", 384/N, logarithmic_scalar_color_map)
	}

	{
		write_scalar(pressure_ping, "pressure.png", 384/N, logarithmic_scalar_color_map)

		residual := make_2D(f32, N, N)
		defer delete_2D(residual)	
		calc_residual(pressure_ping, initial_divergence, residual)	
		write_scalar(residual, "residual.png", 384/N, logarithmic_scalar_color_map)
		write_scalar(residual, "residual_linear.png", 384/N, linear_scalar_color_map, 2.0)

		velocity := make_2D([2]f32, N, N)
		defer delete_2D(velocity)
		subtract_gradient(pressure_ping, initial_velocity, velocity)
		write_velocity(velocity, "velocity.png", 384/N)

		divergence := make_2D(f32, N, N)
		defer delete_2D(divergence)
		calc_divergence(velocity, divergence)
		write_scalar(divergence, "divergence.png", 384/N, logarithmic_scalar_color_map)
		write_scalar(divergence, "divergence_linear.png", 384/N, linear_scalar_color_map)

		max_residual, rms_residual, gradient_residual := calc_stats(residual)

		fmt.printf("max residual:       %e %e\n", max_residual, max_residual/max_residual0)
		fmt.printf("quadratic residual: %e %e\n", rms_residual, rms_residual/rms_residual0)
		fmt.printf("gradient residual:  %e %e\n", gradient_residual, gradient_residual/gradient_residual0)
	}
}

calc_residual :: proc(x, b, r: [][]f32) {
	N := len(x)
	assert(len(x[0]) == N)
	assert(len(b) == N)
	assert(len(b[0]) == N)
	assert(len(r) == N)
	assert(len(r[0]) == N)

	// We explicitly skip the lower border, which is initialized to 0 already
	// We implicitly skip the upper border, which is implicitly defined to be 0
	for j in 1..<N {
		for i in 1..<N {
			center      := x[j+0][i+0]
			lower_left  := x[j-1][i-1]
			lower_right := x[j-1][i+1] if i+1 < N else 0.0
			upper_left  := x[j+1][i-1] if j+1 < N else 0.0
			upper_right := x[j+1][i+1] if i+1 < N && j+1 < N else 0.0

			r[j][i] = b[j][i] + lower_left + lower_right + upper_left + upper_right - 4*center
		}
	}
}

iterate :: proc(in_x, b, out_x: [][]f32, omega: f32 = 1.0) {
	N := len(in_x)
	assert(len(in_x[0]) == N)
	assert(len(b) == N)
	assert(len(b[0]) == N)
	assert(len(out_x) == N)
	assert(len(out_x[0]) == N)

	// We explicitly skip the lower border, which is initialized to 0 already
	// We implicitly skip the upper border, which is implicitly defined to be 0
	for j in 1..<N {
		for i in 1..<N {
			lower_left  := in_x[j-1][i-1]
			lower_right := in_x[j-1][i+1] if i+1 < N else 0.0
			upper_left  := in_x[j+1][i-1] if j+1 < N else 0.0
			upper_right := in_x[j+1][i+1] if i+1 < N && j+1 < N else 0.0

			// 2.0 = sqrt(2) * sqrt(2) is the diagonal grid spacing
			x := (b[j][i] + lower_left + lower_right + upper_left + upper_right) / 4.0				

			// Weighted Jacobi if in_x != out_x and omega <= 1
			// SOR if in_x == out_x and 1 <= omega <= 2
			out_x[j][i] = in_x[j][i] * (1.0 - omega) + x * omega
		}
	}
}

subtract_gradient :: proc(pressure: [][]f32, in_velocity, out_velocity: [][][2]f32) {
	N := len(pressure)
	assert(len(pressure[0]) == N)
	assert(len(in_velocity) == N)
	assert(len(in_velocity[0]) == N)
	assert(len(out_velocity) == N)
	assert(len(out_velocity[0]) == N)

	for j in 0..<N {
		for i in 0..<N {
			lower_left  := pressure[j+0][i+0]
			lower_right := pressure[j+0][i+1] if i+1 < N else 0.0
			upper_left  := pressure[j+1][i+0] if j+1 < N else 0.0
			upper_right := pressure[j+1][i+1] if i+1 < N && j+1 < N else 0.0

			left  := (lower_left  + upper_left)  / 2.0
			right := (lower_right + upper_right) / 2.0
			up    := (upper_left  + upper_right) / 2.0
			down  := (lower_left  + lower_right) / 2.0

			out_velocity[j][i] = in_velocity[j][i] - {
				right - left,
				up - down,
			}
		}
	}
}

calc_divergence :: proc(velocity: [][][2]f32, divergence: [][]f32) {
	N := len(velocity)
	assert(len(velocity[0]) == N)
	assert(len(divergence) == N)
	assert(len(divergence[0]) == N)

	// We explicitly skip the lower border, which is initialized to 0 already
	// We implicitly skip the upper border, which is implicitly defined to be 0
	for j in 1..<N {
		for i in 1..<N {
			// Fetch the nearest 4 velocities
			lower_left  := velocity[j-1][i-1]
			lower_right := velocity[j-1][i+0]
			upper_left  := velocity[j+0][i-1]
			upper_right := velocity[j+0][i+0]

			// Calculate the face values
			left  := (lower_left.x  + upper_left.x)  / 2.0
			right := (lower_right.x + upper_right.x) / 2.0
			up    := (upper_left.y  + upper_right.y) / 2.0
			down  := (lower_left.y  + lower_right.y) / 2.0

			// Central difference using face values, grid spacing 1
			// 2.0 = sqrt(2) * sqrt(2) is the diagonal grid spacing
			divergence[j][i] = (-math.SQRT_TWO*math.SQRT_TWO)*(right - left + up - down)
		}
	}
}


/*  
	Vertex Grid 2D:
		- Velocity at cell center
		- Divergence/Pressure at cell corners


	Fine Grid: 9x5 Pressure Points, 8x4 Velocity Cells

	 P ───────── P ───────── P ───────── P ───────── P ───────── P ───────── P ───────── P ───────── P
	 ╷           ╷           ╷           ╷           ╷           ╷           ╷           ╷           ╷
	 │           │           │           │           │           │           │           │           │
	 │    x,y    │    x,y    │    x,y    │    x,y    │    x,y    │    x,y    │    x,y    │    x,y    │
	 │           │           │           │           │           │           │           │           │
	 ╵           ╵           ╵           ╵           ╵           ╵           ╵           ╵           ╵
	 P ───────── p ───────── p ───────── p ───────── p ───────── p ───────── p ───────── p ───────── P
	 ╷           ╷           ╷           ╷           ╷           ╷           ╷           ╷           ╷
	 │           │           │           │           │           │           │           │           │
	 │    x,y    │    x,y    │    x,y    │    x,y    │    x,y    │    x,y    │    x,y    │    x,y    │
	 │           │           │           │           │           │           │           │           │
	 ╵           ╵           ╵           ╵           ╵           ╵           ╵           ╵           ╵
	 P ───────── p ───────── p ───────── p ───────── p ───────── p ───────── p ───────── p ───────── P
	 ╷           ╷           ╷           ╷           ╷           ╷           ╷           ╷           ╷
	 │           │           │           │           │           │           │           │           │
	 │    x,y    │    x,y    │    x,y    │    x,y    │    x,y    │    x,y    │    x,y    │    x,y    │
	 │           │           │           │           │           │           │           │           │
	 ╵           ╵           ╵           ╵           ╵           ╵           ╵           ╵           ╵
	 P ───────── p ───────── p ───────── p ───────── p ───────── p ───────── p ───────── p ───────── P
	 ╷           ╷           ╷           ╷           ╷           ╷           ╷           ╷           ╷
	 │           │           │           │           │           │           │           │           │
	 │    x,y    │    x,y    │    x,y    │    x,y    │    x,y    │    x,y    │    x,y    │    x,y    │
	 │           │           │           │           │           │           │           │           │
	 ╵           ╵           ╵           ╵           ╵           ╵           ╵           ╵           ╵
	 P ───────── p ───────── p ───────── p ───────── p ───────── p ───────── p ───────── p ───────── P

	
	Coarse Grid: 5x3 Pressure Points

	 P ───────────────────── P ───────────────────── P ───────────────────── P ───────────────────── P
	 ╷                       ╷                       ╷                       ╷                       ╷
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 ╵                       ╵                       ╵                       ╵                       ╵
	 P ───────────────────── P ───────────────────── P ───────────────────── P ───────────────────── P
	 ╷                       ╷                       ╷                       ╷                       ╷
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 ╵                       ╵                       ╵                       ╵                       ╵
	 P ───────────────────── P ───────────────────── P ───────────────────── P ───────────────────── P



	Indexing:

	0,4 ─────── 1,4 ─────── 2,4 ─────── 3,4 ─────── 4,4 ─────── 5,4 ─────── 6,4 ─────── 7,4 ─────── 8,4
	 ╷           ╷           ╷           ╷           ╷           ╷           ╷           ╷           ╷
	 │           │           │           │           │           │           │           │           │
	 │    0,3    │    1,3    │    2,3    │    3,3    │    4,3    │    5,3    │    6,3    │    7,3    │
	 │           │           │           │           │           │           │           │           │
	 ╵           ╵           ╵           ╵           ╵           ╵           ╵           ╵           ╵
	0,3 ─────── 1,3 ─────── 2,3 ─────── 3,3 ─────── 4,3 ─────── 5,3 ─────── 6,3 ─────── 7,3 ─────── 8,3
	 ╷           ╷           ╷           ╷           ╷           ╷           ╷           ╷           ╷
	 │           │           │           │           │           │           │           │           │
	 │    0,2    │    1,2    │    2,2    │    3,2    │    4,2    │    5,2    │    6,2    │    7,2    │
	 │           │           │           │           │           │           │           │           │
	 ╵           ╵           ╵           ╵           ╵           ╵           ╵           ╵           ╵
	0,2 ─────── 1,2 ─────── 2,2 ─────── 3,2 ─────── 4,2 ─────── 5,2 ─────── 6,2 ─────── 7,2 ─────── 8,2
	 ╷           ╷           ╷           ╷           ╷           ╷           ╷           ╷           ╷
	 │           │           │           │           │           │           │           │           │
	 │    0,1    │    1,1    │    2,1    │    3,1    │    4,1    │    5,1    │    6,1    │    7,1    │
	 │           │           │           │           │           │           │           │           │
	 ╵           ╵           ╵           ╵           ╵           ╵           ╵           ╵           ╵
	0,1 ─────── 1,1 ─────── 2,1 ─────── 3,1 ─────── 4,1 ─────── 5,1 ─────── 6,1 ─────── 7,1 ─────── 8,1
	 ╷           ╷           ╷           ╷           ╷           ╷           ╷           ╷           ╷
	 │           │           │           │           │           │           │           │           │
	 │    0,0    │    1,0    │    2,0    │    3,0    │    4,0    │    5,0    │    6,0    │    7,0    │
	 │           │           │           │           │           │           │           │           │
	 ╵           ╵           ╵           ╵           ╵           ╵           ╵           ╵           ╵
	0,0 ─────── 1,0 ─────── 2,0 ─────── 3,0 ─────── 4,0 ─────── 5,0 ─────── 6,0 ─────── 7,0 ─────── 8,0


	Coarsening:

	16*pc[j][i] = 1*(p[2*j-1][2*i-1] + p[2*j-1][2*i+1] + p[2*j+1][2*i-1] + p[2*j+1][2*i+1])
	            + 2*(p[2*j-1][2*i+0] + p[2*j+0][2*i-1] + p[2*j+0][2*i+1] + p[2*j+1][2*i+0])
	            + 4*p[2*j][2*i]

	0,2 ─────────────────── 1,2 ─────────────────── 2,2 ─────────────────── 3,2 ─────────────────── 4,2
	 ╷                       ╷                       ╷                       ╷                       ╷
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 ╵                       ╵                       ╵                       ╵                       ╵
	0,1 ─────────────────── 1,1 ─────────────────── 2,1 ─────────────────── 3,1 ─────────────────── 4,1
	 ╷                       ╷                       ╷                       ╷                       ╷
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 │                       │                       │                       │                       │
	 ╵                       ╵                       ╵                       ╵                       ╵
	0,0 ─────────────────── 1,0 ─────────────────── 2,0 ─────────────────── 3,0 ─────────────────── 4,0


	Interpolating:

	x := f32(i % 2) / 2.0
	y := f32(j % 2) / 2.0
	pc00 := pc[j/2+0][i/2+0]
	pc10 := pc[j/2+0][i/2+1] if i/2 == M/2 else 0.0
	pc01 := pc[j/2+1][i/2+0] if j/2 == N/2 else 0.0
	pc11 := pc[j/2+1][i/2+1] if i/2 == M/2 || j/2 == N/2 else 0.0
	p[j][i] = pc00 * (1 - x) * (1 - y)
	        + pc10 * x * (1 - y)
	        + pc01 * (1 - x) * y
	        + pc11 * x * y

*/

restrict :: proc(in_r, out_r: [][]f32) {
	N := len(out_r)
	assert(N == len(out_r[0]))
	assert(2*N == len(in_r))
	assert(2*N == len(in_r[0]))

	for j in 1..<N {
		for i in 1..<N {
			corners := in_r[2*j-1][2*i-1] + in_r[2*j-1][2*i+1] + in_r[2*j+1][2*i-1] + in_r[2*j+1][2*i+1]
			edges   := in_r[2*j-1][2*i+0] + in_r[2*j+0][2*i-1] + in_r[2*j+0][2*i+1] + in_r[2*j+1][2*i+0]
			center  := in_r[2*j+0][2*i+0]

			r := (4.0 * center + 2.0 * edges + 1.0 * corners) / 16.0
			//r = (4.0 * center + 4.0 * corners) / 20.0
			//r = (4.0 * center + 1.0 * corners) / 8.0
			//r = center;
			out_r[j][i] = 4.0 * r
		}
	}
}

prolongate :: proc(in_error, out_error: [][]f32) {
	N := len(out_error)
	assert(N == len(out_error[0]))
	assert(N == 2*len(in_error))
	assert(N == 2*len(in_error[0]))

	for j in 1..<N {
		for i in 1..<N {
			x := f32(i % 2) / 2.0
			y := f32(j % 2) / 2.0
			e00 := in_error[j/2+0][i/2+0]
			e10 := in_error[j/2+0][i/2+1] if i/2+1 != N/2 else 0.0
			e01 := in_error[j/2+1][i/2+0] if j/2+1 != N/2 else 0.0
			e11 := in_error[j/2+1][i/2+1] if i/2+1 != N/2 && j/2+1 != N/2 else 0.0

			out_error[j][i] = e00 * (1 - x) * (1 - y) \
			                + e10 * x * (1 - y) \
			                + e01 * (1 - x) * y \
			                + e11 * x * y
		}
	}
}



vcycle :: proc(x_ping, rhs, x_pong: ^[][]f32, pre_smooths, solves, post_smooths: int) {
	N := len(x_ping)
	assert(N == len(x_ping[0]))
	assert(N == len(x_pong))
	assert(N == len(x_pong[0]))
	assert(N == len(rhs))
	assert(N == len(rhs[0]))
	//fmt.println("N =", N)
	if N/2 <= 1 {
	//if N == MAIN_N/2 {
		//fmt.println("Solve")
		// solve
		for iter in 0..<4 {
			iterate(x_ping^, rhs^, x_ping^, 1.0) // Gauss-Seidel
		}
		write_scalar(x_ping^, "pressure_coarse.png", 384/N, logarithmic_scalar_color_map)

		residual := x_pong
		calc_residual(x_ping^, rhs^, residual^)
		write_scalar(residual^, "residual_coarse2.png", 384/N, logarithmic_scalar_color_map)
		return
	}
	//fmt.println("Deeper")

	for iter := 0; iter < 2; iter += 1 {
		iterate(x_ping^, rhs^, x_pong^, 8/10.0) // Weighted Jacobi
		x_ping^, x_pong^ = x_pong^, x_ping^
	}
	residual := x_pong
	calc_residual(x_ping^, rhs^, residual^)
	write_scalar(residual^, "residual_fine.png", 384/N, logarithmic_scalar_color_map)
	write_scalar(residual^, "residual_fine_linear.png", 384/N, linear_scalar_color_map, 2.0)

	residual_coarse := make_2D(f32, N/2, N/2, context.temp_allocator)
	restrict(residual^, residual_coarse)
	write_scalar(residual_coarse, "residual_coarse.png", 2*384/N, logarithmic_scalar_color_map)
	write_scalar(residual_coarse, "residual_coarse_linear.png", 2*384/N, linear_scalar_color_map, 2.0)
	
	x_ping_coarse := make_2D(f32, N/2, N/2, context.temp_allocator)
	x_pong_coarse := make_2D(f32, N/2, N/2, context.temp_allocator)
	vcycle(&x_ping_coarse, &residual_coarse, &x_pong_coarse, pre_smooths, solves, post_smooths)

	e_fine := x_pong
	prolongate(x_ping_coarse, e_fine^)
	write_scalar(x_ping_coarse, "error_coarse.png", 2*384/N, logarithmic_scalar_color_map)
	write_scalar(e_fine^, "error_fine.png", 384/N, logarithmic_scalar_color_map)
	for j in 1..<N do for i in 1..<N do x_ping[j][i] = x_ping[j][i] + e_fine[j][i]
	write_scalar(x_ping^, "pressure_fine.png", 384/N, logarithmic_scalar_color_map)

	//for j in 0..<N do fmt.printf("%+f\n", x_ping[j])


	//for iter := 0; iter < post_smooths; iter += 1 {
	calc_residual(x_ping^, rhs^, x_pong^)
	max_residual0, rms_residual0, gradient_residual0 := calc_stats(x_pong^)
	//fmt.printf("   residual: %.6f %.6f %.6f\n", max_residual0, rms_residual0, gradient_residual0)
	//fmt.printf("%d: residual: %.6f %.6f %.6f\n", 0, max_residual0/max_residual0, rms_residual0/rms_residual0, gradient_residual0/gradient_residual0)

	for iter := 0; iter < 4; iter += 1 {
		iterate(x_ping^, rhs^, x_pong^, 8.0/10.0) // Weighted Jacobi
		x_ping^, x_pong^ = x_pong^, x_ping^

		calc_residual(x_ping^, rhs^, x_pong^)
		max_residual, rms_residual, gradient_residual := calc_stats(x_pong^)
		//fmt.printf("%d: residual: %.6f %.6f %.6f\n", iter+1, max_residual/max_residual0, rms_residual/rms_residual0, gradient_residual/gradient_residual0)
	}
}

ITER_MAX :: 3


calc_stats :: proc(grid: [][]f32) -> (max_value, rms_value, gradient_value: f32) {
	N := len(grid)
	assert(N == len(grid[0]))

	for j in 1..<N {
		for i in 1..<N {
			max_value = max(max_value, abs(grid[j][i]))
			rms_value += grid[j][i] * grid[j][i]
		}
	}

	for j in 2..<N {
		for i in 2..<N {
			dx := grid[j][i] - grid[j][i-1]
			dy := grid[j][i] - grid[j-1][i]
			gradient_value += dx*dx + dy*dy
		}
	}
	rms_value = math.sqrt(rms_value / f32((N-1) * (N-1)))
	gradient_value = math.sqrt(gradient_value / f32((N-2) * (N-2)))
	return
} 
