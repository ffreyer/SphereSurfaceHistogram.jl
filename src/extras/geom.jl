"""
    to_cartesian(theta, phi)

Calculates the cartesian vector (x, y, z) for a given pair of angles (polar/z
angle theta, azimut/xy angle phi).
"""
function to_cartesian(theta, phi)
    Vec3f0(
        sin(theta) * cos(phi),
        sin(theta) * sin(phi),
        cos(theta)
    )
end


"""
    face_mesh(B::SSHBinner, extrude=0)

Generates a tiled/disconnected mesh where each face represents a bin.

`extrude` allows faces to be moved outwards (≈ (1 + extrude)eᵣ). For
`extrude > 0` the faces will sperate, making their shape easier to see.
"""
function face_mesh(B::SSHBinner, extrude=0)
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
    vertex_mesh(B::SSHBinner)

Generates a mesh where each vertex represents (the center of) a bin.

This type of mesh makes it easy to visualize the histogram with colors. Mapping
each bin-filling to a color (e.g. using `to_hue(B)`) gives colors
corresponding to the vertices of this mesh (in order). This can be plotted with
`Makie` using `Makie.mesh(vertexl_mesh(B), color = colors)`.

See also: [`bin_positions`](@ref), [`face_mesh`](@ref)
"""
function vertex_mesh(B::SSHBinner)
    phis = [
        #   lvl-jump     lvl-step
        2pi * (i-1 + (j-0.5)/B.phi_divisions[i])
        for i in eachindex(B.phi_divisions)
            for j in 1:B.phi_divisions[i]
    ]
    phis[1] += 0.5pi # increasing the initial angle fixes overlapping triangles
    phis[end] += 0.5pi

    # 2N - 4
    faces = Face{3, Int}[]
    j = 1
    N = length(phis)
    for i in 3:N
        push!(faces, Face{3, Int}(j, i-1, i))
        while phis[i] > phis[j] + 2pi
            push!(faces, Face{3, Int}(j, i, j+1))
            j += 1
        end
    end
    push!(faces, Face{3, Int}(N-2, N-1, N))

    GLNormalMesh(vertices = dual_points(B), faces = faces)
end


"""
    bin_positions(B::SSHBinner)

Returns a list of points corresponding to the center positions of each bin.
This list is equivalent to the vertices used in `vertex_mesh(B)`.

See also: [`vertex_mesh`](@ref)
"""
function bin_positions(B::SSHBinner)
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


"""
    line_segments_minimal(B::SSHBinner)

Generates linesegments marking each bin area. This version is minimal in the
sense that doesn't interpolate any points.

See also: [`line_segments`](@ref)
"""
function line_segments_minimal(B::SSHBinner)
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


"""
    line_segments(B::SSHBinner; N_fragments = 32)

Generates line segments marking each bin area. Extra points may be interpolated
to generate smooth arcs.

Specifically, points are added to an arc, such that it would include at least
`N_fragments` points over an angle 0:pi.

See also: [`line_segments_minimal`](@ref)
"""
function line_segments(B::SSHBinner; N_fragments=32)
    lines = Point3f0[]

    # horizontal rings
    for theta in B.thetas[2:end-1]
        z = cos(theta)
        # Number of steps
        R_xy = sqrt(1.0 - z^2)
        # circumference = 2πR, make more points if longer line
        for k in 0:2N_fragments-1
            phi1 = pi * k / N_fragments
            phi2 = pi * (k+1) / N_fragments
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
