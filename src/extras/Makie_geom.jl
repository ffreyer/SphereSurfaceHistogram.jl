"""
    to_rects(B::SSHBinner, extrude=0)

Generates a mesh where each face represents a bin.

`extrude` allows faces to be moved outwards (≈ (1 + extrude)eᵣ). For
`extrude > 0` the faces will sperate, makign their shape easier to see.
"""
function to_rects(B::SSHBinner, extrude=0)
    vertices = Point3f0[]
    faces = Face{3, Int}[]

    push!(vertices, Point3f0(0, 0, 1) + extrude * Point3f0(0,0,1))
    theta2 = B.thetas[2]
    N = B.phi_divisions[2]
    for k in 0:N-1
        phi1 = 2pi * k/N
        phi2 = 2pi * (k+1) / N
        pos_1 = to_cartesian(theta2, phi1) + extrude * Point3f0(0,0,1)
        pos_2 = to_cartesian(theta2, phi2) + extrude * Point3f0(0,0,1)
        push!(vertices, pos_1)
        push!(vertices, pos_2)
        push!(faces, Face(1, 2+2k, 3+2k))
    end

    for i in 2:length(B.phi_divisions)-1
        theta1 = B.thetas[i]
        theta2 = B.thetas[i+1]
        N = B.phi_divisions[i]

        for k in 0:N-1
            phi1 = 2pi * k / N
            phi2 = 2pi * (k+1) / N
            pos_11 = to_cartesian(theta1, phi1)
            pos_12 = to_cartesian(theta1, phi2)
            pos_21 = to_cartesian(theta2, phi1)
            pos_22 = to_cartesian(theta2, phi2)

            # extrude "rects"
            dir = 0.25(pos_11 + pos_12 + pos_21 + pos_22)
            pos_11 += extrude*dir
            pos_12 += extrude*dir
            pos_21 += extrude*dir
            pos_22 += extrude*dir

            l = length(vertices)
            append!(vertices, [pos_11, pos_12, pos_21, pos_22])
            push!(faces, Face(l+1, l+3, l+2))
            push!(faces, Face(l+2, l+3, l+4))
        end
    end

    push!(vertices, Point3f0(0, 0, -1) + extrude * Point3f0(0,0,-1))
    l = length(vertices)
    theta2 = B.thetas[end-1]
    N = B.phi_divisions[end-1]

    for k in 0:N-1
        phi1 = 2pi * k / N
        phi2 = 2pi * (k+1) / N
        pos_1 = to_cartesian(theta2, phi1) + extrude * Point3f0(0,0,-1)
        pos_2 = to_cartesian(theta2, phi2) + extrude * Point3f0(0,0,-1)
        push!(vertices, pos_1)
        push!(vertices, pos_2)
        push!(faces, Face(l, l+2+2k, l+1+2k))
    end

    GLNormalMesh(vertices = vertices, faces = faces)
end


"""
    to_dual_mesh(B::SSHBinner)

Generates a mesh where each vertex represents (the center of) a bin.

This type of mesh makes it easy to visualize the histogram with colors. Mapping
each bin-filling to a color (e.g. using `to_hue(B)`) gives colors
corresponding to the vertices of this mesh (in order). This can be plotted with
`Makie` using `Makie.mesh(to_dual_mesh(B), color = colors)`.
"""
function to_dual_mesh(B::SSHBinner)
    vertices = Point3f0[]
    faces = Face{3, Int}[]

    # thetas[1]:thetas[2] are 1 area = phi_divs[1]          center @ 0, 0, 1
    # thetas[2]:thetas[3] are N areas = phi_divs[2]
    # thetas[l-2]:thetas[l-1] are N' areas = phi_divs[l-2]
    # thetas[l-1]:thetas[l] are 1 area = phi_divs[l-1]      center @ 0, 0, -1


    # This is the center between thetas[1] and thetas[2]!
    # The representation becomes weird though....
    push!(vertices, Point3f0(0, 0, 1))

    # connecting 1 to N doesn't work with the method used for N to M
    # same with N to 1, see below
    theta = 0.5(B.thetas[2] + B.thetas[3])
    N = B.phi_divisions[2]
    for k in 0:N-1
        phi = 2pi * (k + 0.5) / N
        push!(vertices, to_cartesian(theta, phi))
        push!(faces, Face(1, k+2, ((k+1)%N)+2))
    end

    for i in 3:length(B.thetas)-2
        theta = 0.5(B.thetas[i] + B.thetas[i+1])
        N = B.phi_divisions[i]

        for k in 0:N-1
            phi = 2pi * (k + 0.5) / N
            push!(vertices, to_cartesian(theta, phi))

            l = length(vertices)
            l < 3 && continue
            # I should probably find a better name that M and N
            M = B.phi_divisions[i-1]

            if M == N
                push!(faces, Face(l-1, l, l-N))
                push!(faces, Face(l, l-N+1, l-N))

            elseif M < N
                # offset to top layer
                # if factor = 2: increment inner_index every 2nd step
                # if factor = 4: increment inner_index every 4th step
                factor = div(N, M)
                # equivalent to l - phi_divs - div(factor-1 + k*(factor-1), factor)
                inner_index = l - M - div((k+1)*(factor-1), factor)
                # println("l = $l   \t  inner_index = $inner_index")
                if N-1 > k > 0
                    if k % 2 == 1
                        push!(faces, Face(l-1, l, inner_index))
                        push!(faces, Face(l, inner_index+1, inner_index))
                    else
                        push!(faces, Face(l-1, l, inner_index))
                    end
                elseif k == 0
                    last_index = l + N - 1
                    push!(faces, Face(last_index-1, last_index, l-1))
                    push!(faces, Face(last_index, inner_index, l-1))
                    push!(faces, Face(last_index, l, inner_index))
                end

            elseif M > N
                offset = M - k
                if k > 0
                    push!(faces, Face(l-1, l, l-offset-1))
                    push!(faces, Face(l, l-offset, l-offset-1))
                    push!(faces, Face(l, l-offset+1, l-offset))
                else
                    push!(faces, Face(l-1, l, l-offset))
                    push!(faces, Face(l, l-offset, l-offset+1))
                    push!(faces, Face(l, l-1, l+N-1))
                end

            end
        end
    end

    # Last layer in loop:
    # thetas[l-2]:thetas[l-1] are N' areas = phi_divs[l-2]  <- this one
    # thetas[l-1]:thetas[l] are 1 area = phi_divs[l-1]      center @ 0, 0, -1
    # so add the last Point and connect it
    push!(vertices, Point3f0(0, 0, -1))
    l = length(vertices)
    N = B.phi_divisions[end-1]
    for k in 0:N-1
        push!(faces, Face(l, l-k-1, l - ((k+1) % N) - 1))
    end

    GLNormalMesh(vertices = vertices, faces = faces)
end


"""
    dual_points(B::SSHBinner)

Generates the same points (vertices) as `to_dual_mesh(B)`, but doesn't generate
the full mesh.
"""
function dual_points(B::SSHBinner)
    points = Point3f0[Point3f0(0, 0, 1)]

    for i in 2:length(B.thetas)-2
        theta = 0.5(B.thetas[i] + B.thetas[i+1])
        N = B.phi_divisions[i]

        for k in 0:N-1
            phi = 2pi * (k + 0.5) / N
            push!(points, to_cartesian(theta, phi))
        end
    end
    push!(points, Point3f0(0, 0, -1))

    points
end


function line_segments(B::SSHBinner)
    lines = Point3f0[]

    # first ring
    theta = B.thetas[2]
    N = B.phi_divisions[2]
    for k in 0:N-1
        phi1 = 2pi * k / N
        phi2 = 2pi * (k+1) / N
        append!(lines, [to_cartesian(theta, phi1), to_cartesian(theta, phi2)])
    end

    # Any thing after that is 1 vertical and 1 horizontal line
    for i in 3:length(B.thetas)-1
        theta1 = B.thetas[i-1]
        theta2 = B.thetas[i]
        if B.phi_divisions[i-1] == B.phi_divisions[i]
            N = B.phi_divisions[i]
        elseif B.phi_divisions[i-1] < B.phi_divisions[i]
            N = B.phi_divisions[i-1]
        else
            N = B.phi_divisions[i-1]
        end

        for k in 0:N-1
            phi1 = 2pi * k / N
            phi2 = 2pi * (k+1) / N

            # horizontal lines
            append!(lines, [to_cartesian(theta2, phi1), to_cartesian(theta2, phi2)])
            # vertical line
            append!(lines, [to_cartesian(theta2, phi1), to_cartesian(theta1, phi1)])
        end
    end

    lines
end


function line_segments2(B::SSHBinner; N_fragments=64)
    lines = Point3f0[]

    # horizontal rings
    for theta in B.thetas[2:end-1]
        z = cos(theta)
        # Number of steps
        R_xy = sqrt(1.0 - z^2)
        # circumference = 2πR, make more points if longer line
        for k in 0:N_fragments-1
            phi1 = 2pi * k / N_fragments
            phi2 = 2pi * (k+1) / N_fragments
            append!(lines, [to_cartesian(theta, phi1), to_cartesian(theta, phi2)])
        end
    end


    # TODO recast this from phi(theta) to theta(phi)
    substeps = div(N_fragments, length(B.thetas))+1
    # vertical lines
    for i in 3:length(B.thetas)-1
        theta1 = B.thetas[i-1]
        theta2 = B.thetas[i]
        if B.phi_divisions[i-1] == B.phi_divisions[i]
            N = B.phi_divisions[i]
        elseif B.phi_divisions[i-1] < B.phi_divisions[i]
            N = B.phi_divisions[i-1]
        else
            N = B.phi_divisions[i-1]
        end

        for k in 0:N-1
            phi1 = 2pi * k / N
            phi2 = 2pi * (k+1) / N
            for c in 0:substeps-1
                t1 = theta1 + c/substeps * (theta2 - theta1)
                t2 = theta1 + (c+1)/substeps * (theta2 - theta1)

                # vertical line
                append!(lines, [
                    to_cartesian(t2, phi1),
                    to_cartesian(t1, phi1)
                ])
            end
        end
    end

    lines
end
