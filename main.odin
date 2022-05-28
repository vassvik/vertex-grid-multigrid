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

main :: proc() {
	N := 32
	assert(N & (N-1) == 0, "N not power of two")

	velocity := make_2D([2]f32, N, N)
	for j in 3*N/8..<5*N/8 {
		fmt.println(j)
		for i in 3*N/8..<5*N/8 {
			velocity[j][i].y = +1.0
		}
	}
	write_velocity(velocity, "velocity.png", 512/N)
}
