package main 

import "core:fmt"
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

write_scalar_linear :: proc(scalar: [][]f32, filename: cstring, SS: int) {
	N := len(scalar)
	assert(len(scalar[0]) == N)
	scalar_image := make_2D([4]u8, SS*N, SS*N)
	defer delete_2D(scalar_image)

	for j in 0..<N {
		for i in 0..<N {
			pixel_value := [4]u8{
				// Flowmap color map
				u8(clamp(-255*scalar[j][i], 0, 255)), // negative values are red
				u8(clamp(+255*scalar[j][i], 0, 255)), // positive values are green
				0,
				255,
			}

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
	iterations := 10
	N := 32
	assert(N & (N-1) == 0, "N not power of two")

	// Initialize Velocity
	velocity := make_2D([2]f32, N, N)
	for j in 3*N/8..<5*N/8 {
		for i in 3*N/8..<5*N/8 {
			velocity[j][i].y = +1.0
		}
	}
	write_velocity(velocity, "initial_velocity.png", 512/N)

	// The pressure and divergence grids are (N+1)x(N+1) in the vertex grid with 
	// NxN velocity cells. We still store these in NxN arrays, where we explicitly
	// maintain the lower boundary value and implicitly assume the upper boundary.
	divergence := make_2D(f32, N, N)
	pressure_ping := make_2D(f32, N, N)
	pressure_pong := make_2D(f32, N, N)

	// Calculate divergence
	// We explicitly skip the lower border, which is initialized to 0 already
	// We implicitly skip the upper border, which is implicitly defined to be 0
	for j in 1..<N {
		for i in 1..<N {
			// Fetch the nearest 4 velocities
			lower_left  := velocity[j-1][i-1]
			lower_right := velocity[j-1][i+0]
			upper_left  := velocity[j+0][i-1]
			upper_right  := velocity[j+0][i+0]

			// Calculate the face values
			left  := (lower_left.x  + upper_left.x) / 2.0
			right := (lower_right.x + upper_right.x) / 2.0
			up    := (upper_left.y  + upper_right.y) / 2.0
			down  := (lower_left.y  + lower_right.y) / 2.0

			// Central difference using face values, grid spacing 1/2
			divergence[j][i] = (right - left) + (up - down)
		}
	}
	write_scalar_linear(divergence, "initial_divergence.png", 512/N)

	// Solve pressure
	// We explicitly skip the lower border, which is initialized to 0 already
	// We implicitly skip the upper border, which is implicitly defined to be 0
	for iter in 0..iterations {
		for j in 1..<N {
			for i in 1..<N {
				lower_left  := pressure_ping[j-1][i-1]
				lower_right := pressure_ping[j-1][i+1] if i+1 < N else 0.0
				upper_left  := pressure_ping[j+1][i-1] if j+1 < N else 0.0
				upper_right := pressure_ping[j+1][i+1] if i+1 < N && j+1 < N else 0.0

				pressure_pong[j][i] = (-2.0*divergence[j][i] + lower_left + lower_right + upper_left + upper_right) / 4.0
			}
		}

		pressure_ping, pressure_pong = pressure_pong, pressure_ping
	}
	write_scalar_linear(pressure_ping, "pressure.png", 512/N)

	// Apply gradient
	for j in 0..<N {
		for i in 0..<N {
			lower_left  := pressure_ping[j+0][i+0]
			lower_right := pressure_ping[j+0][i+1] if i+1 < N else 0.0
			upper_left  := pressure_ping[j+1][i+0] if j+1 < N else 0.0
			upper_right := pressure_ping[j+1][i+1] if i+1 < N && j+1 < N else 0.0

			left  := (lower_left  + upper_left) / 2.0
			right := (lower_right + upper_right) / 2.0
			up    := (upper_left  + upper_right) / 2.0
			down  := (lower_left  + lower_right) / 2.0

			velocity[j][i] -= {
				right - left,
				up - down,
			}
		}

		pressure_ping, pressure_pong = pressure_pong, pressure_ping
	}
	write_velocity(velocity, "velocity.png", 512/N)

	for j in 1..<N {
		for i in 1..<N {
			// Fetch the nearest 4 velocities
			lower_left  := velocity[j-1][i-1]
			lower_right := velocity[j-1][i+0]
			upper_left  := velocity[j+0][i-1]
			upper_right  := velocity[j+0][i+0]

			// Calculate the face values
			left  := (lower_left.x  + upper_left.x) / 2.0
			right := (lower_right.x + upper_right.x) / 2.0
			up    := (upper_left.y  + upper_right.y) / 2.0
			down  := (lower_left.y  + lower_right.y) / 2.0

			// Central difference using face values, grid spacing 1/2
			divergence[j][i] = (right - left) + (up - down)
		}
	}
	write_scalar_linear(divergence, "divergence.png", 512/N)
}
