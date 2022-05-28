package main 

import "core:fmt"
import "core:math"
import stbi "vendor:stb/image"

make_2D :: proc($T: typeid, M, N: int) -> [][]T {
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
    i := math.floor(4 * ((logv + 7.0) - f));

    c: [4]f32 = {0, 0, 0, 1}
         if f < 0.0 do c.rgb = {0.0, 0.0, 0.0};
    else if f < 1.0 do c.rgb = mix({1.0, 0.0, 0.0}, {1.0, 1.0, 1.0}, i / 4.0);
    else if f < 2.0 do c.rgb = mix({0.0, 1.0, 0.0}, {1.0, 1.0, 1.0}, i / 4.0);
    else if f < 3.0 do c.rgb = mix({0.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, i / 4.0);
    else if f < 4.0 do c.rgb = mix({1.0, 1.0, 0.0}, {1.0, 1.0, 1.0}, i / 4.0);
    else if f < 5.0 do c.rgb = mix({1.0, 0.0, 1.0}, {1.0, 1.0, 1.0}, i / 4.0);
    else if f < 6.0 do c.rgb = mix({0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, i / 4.0);
    else if f < 7.0 do c.rgb = mix({1.0, 0.5, 0.0}, {1.0, 1.0, 1.0}, i / 4.0);
    else if f < 8.0 do c.rgb = mix({1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, i / 4.0);

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

write_scalar :: proc(scalar: [][]f32, filename: cstring, SS: int, color_map: proc(f32) -> [4]u8) {
	N := len(scalar)
	assert(len(scalar[0]) == N)
	scalar_image := make_2D([4]u8, SS*N, SS*N)
	defer delete_2D(scalar_image)

	for j in 0..<N {
		for i in 0..<N {
			pixel_value := color_map(scalar[j][i])

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


main :: proc() {
	stbi.flip_vertically_on_write(true)

	N := 32
	dump_interval := N
	iterations := N*N
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

	// Solve pressure
	for iter := 0; iter < iterations; iter += 1 {
		fmt.println(iter)
		iterate_jacobi(pressure_ping, initial_divergence, pressure_pong)
		pressure_ping, pressure_pong = pressure_pong, pressure_ping

		if !((iter % dump_interval) == 0  || iter == iterations-1) do continue
		write_scalar(pressure_ping, "pressure.png", 384/N, logarithmic_scalar_color_map)

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
}

calc_residual :: proc(pressure, divergence, residual: [][]f32) {
	N := len(pressure)
	assert(len(pressure[0]) == N)
	assert(len(divergence) == N)
	assert(len(divergence[0]) == N)
	assert(len(residual) == N)
	assert(len(residual[0]) == N)

	// We explicitly skip the lower border, which is initialized to 0 already
	// We implicitly skip the upper border, which is implicitly defined to be 0
	for j in 1..<N {
		for i in 1..<N {
			center      := pressure[j+0][i+0]
			lower_left  := pressure[j-1][i-1]
			lower_right := pressure[j-1][i+1] if i+1 < N else 0.0
			upper_left  := pressure[j+1][i-1] if j+1 < N else 0.0
			upper_right := pressure[j+1][i+1] if i+1 < N && j+1 < N else 0.0

			// 2.0 = sqrt(2) * sqrt(2) is the diagonal grid spacing
			residual[j][i] = -2.0*divergence[j][i] + lower_left + lower_right + upper_left + upper_right - 4*center
		}
	}
}

iterate_jacobi :: proc(in_pressure, divergence, out_pressure: [][]f32, omega: f32 = 1.0) {
	N := len(in_pressure)
	assert(len(in_pressure[0]) == N)
	assert(len(divergence) == N)
	assert(len(divergence[0]) == N)
	assert(len(out_pressure) == N)
	assert(len(out_pressure[0]) == N)

	// We explicitly skip the lower border, which is initialized to 0 already
	// We implicitly skip the upper border, which is implicitly defined to be 0
	for j in 1..<N {
		for i in 1..<N {
			lower_left  := in_pressure[j-1][i-1]
			lower_right := in_pressure[j-1][i+1] if i+1 < N else 0.0
			upper_left  := in_pressure[j+1][i-1] if j+1 < N else 0.0
			upper_right := in_pressure[j+1][i+1] if i+1 < N && j+1 < N else 0.0

			// 2.0 = sqrt(2) * sqrt(2) is the diagonal grid spacing
			p_Jacobi := (-2.0*divergence[j][i] + lower_left + lower_right + upper_left + upper_right) / 4.0				

			// Weighted Jacobi if in_pressure != out_pressure and omega <= 1
			// SOR if in_pressure == out_pressure and 1 <= omega <= 2
			out_pressure[j][i] = in_pressure[j][i] * (1.0 - omega) + p_Jacobi * omega
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

			// Central difference using face values, grid spacing 1/2
			divergence[j][i] = (right - left) + (up - down)
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
	p[j][i] = pc[j/2+0][i/2+0] * (1 - x) * (1 - y)
	        + pc[j/2+1][i/2+0] * x * (1 - y)
	        + pc[j/2+0][i/2+1] * (1 - x) * y
	        + pc[j/2+1][i/2+1] * x * y

*/