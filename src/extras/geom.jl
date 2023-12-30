"""
    to_cartesian(theta, phi)

Calculates the cartesian vector (x, y, z) for a given pair of angles (polar/z
angle theta, azimut/xy angle phi).
"""
function to_cartesian(radius, theta, phi)
    return Vec3f(
        radius * sin(theta) * cos(phi),
        radius * sin(theta) * sin(phi),
        radius * cos(theta)
    )
end
function to_cartesian(theta, phi)
    return Vec3f(
        sin(theta) * cos(phi),
        sin(theta) * sin(phi),
        cos(theta)
    )
end


# TODO: this has holes
"""
    face_mesh(B::AbstractSSH, radius = 1, extrude = 0)

Generates a tiled/disconnected mesh where each face represents a bin.

`extrude` allows faces to be moved outwards (≈ (1 + extrude)eᵣ). For
`extrude > 0` the faces will sperate, making their shape easier to see.
"""
face_mesh(B::AbstractSSH, radius = 1, extrude = 0) = face_mesh(B.tessellation, radius, extrude)
function face_mesh(B::SphereTessellationMap, radius, extrude)
    vertices = Point3f[]
    faces = GLTriangleFace[]

    push!(vertices, Point3f(0, 0, radius) + extrude * Point3f(0,0,1))
    theta2 = B.thetas[2]
    N = B.phi_divisions[2]
    for k in 0:N-1
        phi1 = 2pi * k/N
        phi2 = 2pi * (k+1) / N
        pos_1 = to_cartesian(radius, theta2, phi1) + extrude * Point3f(0,0,1)
        pos_2 = to_cartesian(radius, theta2, phi2) + extrude * Point3f(0,0,1)
        push!(vertices, pos_1)
        push!(vertices, pos_2)
        push!(faces, GLTriangleFace(1, 2+2k, 3+2k))
    end

    for i in 2:length(B.phi_divisions)-1
        theta1 = B.thetas[i]
        theta2 = B.thetas[i+1]
        N = B.phi_divisions[i]

        for k in 0:N-1
            phi1 = 2pi * k / N
            phi2 = 2pi * (k+1) / N
            pos_11 = to_cartesian(radius, theta1, phi1)
            pos_12 = to_cartesian(radius, theta1, phi2)
            pos_21 = to_cartesian(radius, theta2, phi1)
            pos_22 = to_cartesian(radius, theta2, phi2)

            # extrude "rects"
            dir = 0.25(pos_11 + pos_12 + pos_21 + pos_22)
            pos_11 += extrude*dir
            pos_12 += extrude*dir
            pos_21 += extrude*dir
            pos_22 += extrude*dir

            l = length(vertices)
            push!(vertices, pos_11, pos_12, pos_21, pos_22)
            push!(faces, GLTriangleFace(l+1, l+3, l+2))
            push!(faces, GLTriangleFace(l+2, l+3, l+4))
        end
    end

    push!(vertices, Point3f(0, 0, -1) + extrude * Point3f(0,0,-1))
    l = length(vertices)
    theta2 = B.thetas[end-1]
    N = B.phi_divisions[end-1]

    for k in 0:N-1
        phi1 = 2pi * k / N
        phi2 = 2pi * (k+1) / N
        pos_1 = to_cartesian(radius, theta2, phi1) + extrude * Point3f(0,0,-1)
        pos_2 = to_cartesian(radius, theta2, phi2) + extrude * Point3f(0,0,-1)
        push!(vertices, pos_1)
        push!(vertices, pos_2)
        push!(faces, GLTriangleFace(l, l+2+2k, l+1+2k))
    end

    # GeometryBasics.Mesh(GeometryBasics.meta(vertices; normals=ns), faces)
    return normal_mesh(vertices, faces)
end

# TODO: this is not correct - faces missplaced
function voxel_mesh(B::AbstractSSH; radius = 1, thickness = 0.1)
    return voxel_mesh(B.tessellation, radius, thickness)
end
function voxel_mesh(B::SphereTessellationMap, radius, thickness)
    _mesh = face_mesh(B, radius)
    vs = decompose(Point3f, _mesh)
    fs = faces(_mesh)
    N = length(vs)

    voxel_vs = vcat(
        (radius - 0.5 * thickness) / radius .* vs,
        (radius + 0.5 * thickness) / radius .* vs
    )
    voxel_fs = vcat(fs, [GLTriangleFace(f[1] + N, f[2] + N, f[3] + N) for f in fs])

    # top face
    # exclude center point at (0, 0, z)
    N = B.phi_divisions[2]
    for k in 0:N-1
        push!(voxel_fs, GLTriangleFace(2+2k, 3+2k, 2+2k+N), GLTriangleFace(3+2k, 3+2k+N, 2+2k+N))
    end

    # vertex counter
    l = 2*(N-1) + 1

    # main body faces
    # each section here is a quad
    for i in 2:length(B.phi_divisions)-1
        N = B.phi_divisions[i]
        for k in 0:N-1
            l += 4
            push!(voxel_fs,
                GLTriangleFace(l+1, l+3, l+1+N), GLTriangleFace(l+3, l+3+N, l+1+N),
                GLTriangleFace(l+3, l+4, l+3+N), GLTriangleFace(l+4, l+4+N, l+3+N),
                GLTriangleFace(l+4, l+2, l+4+N), GLTriangleFace(l+2, l+2+N, l+4+N),
                GLTriangleFace(l+2, l+1, l+2+N), GLTriangleFace(l+1, l+1+N, l+2+N)
            )
        end
    end

    # l+1 is the center, exclude that
    N = B.phi_divisions[end-1]
    for k in 0:N-1
        push!(voxel_fs, GLTriangleFace(l+3+2k, l+2+2k, l+3+2k+N), GLTriangleFace(l+2+2k, l+2+2k+N, l+3+2k+N))
    end

    return normal_mesh(voxel_vs, voxel_fs)
end


"""
    vertex_mesh(B::AbstractSSH)

Generates a mesh where each vertex represents (the center of) a bin.

This type of mesh makes it easy to visualize the histogram with colors. Mapping
each bin-filling to a color (e.g. using `to_hue(B)`) gives colors
corresponding to the vertices of this mesh (in order). This can be plotted with
`Makie` using `Makie.mesh(vertexl_mesh(B), color = colors)`.

See also: [`bin_positions`](@ref), [`face_mesh`](@ref)
"""
vertex_mesh(B::AbstractSSH, radius = 1.0) = vertex_mesh(B.tessellation, radius)
function vertex_mesh(B::SphereTessellationMap, radius)
    faces = GLTriangleFace[]

    N = 0
    for i in 2:length(B.phi_divisions)
        top_idx = 1
        bot_idx = 1
        top = top_step = 1 / B.phi_divisions[i-1]
        bot = bot_step = 1 / B.phi_divisions[i]
        N2 = N + B.phi_divisions[i-1]

        while (bot < 0.99999) || (top < 0.99999)
            if top < bot
                top += top_step
                top_idx += 1
                i1 = N + mod1(top_idx, B.phi_divisions[i-1])
                i2 = N + top_idx - 1
                i3 = N2 + bot_idx
                push!(faces, GLTriangleFace(i1, i2, i3))
            else
                bot += bot_step
                bot_idx += 1
                i1 = N + top_idx
                i2 = N2 + bot_idx - 1
                i3 = N2 + mod1(bot_idx, B.phi_divisions[i])
                push!(faces, GLTriangleFace(i1, i2, i3))
            end
        end

        if i < length(B.phi_divisions)
            bot_idx += 1
            i1 = N + top_idx
            i2 = N2 + bot_idx - 1
            i3 = N2 + mod1(bot_idx, B.phi_divisions[i])
            push!(faces, GLTriangleFace(i1, i2, i3))
        end

        if i-1 > 1
            top_idx += 1
            i1 = N + mod1(top_idx, B.phi_divisions[i-1])
            i2 = N + top_idx - 1
            i3 = N2 + mod1(bot_idx, B.phi_divisions[i])
            top_idx = mod1(top_idx + 1, B.phi_divisions[i-1])
            push!(faces, GLTriangleFace(i1, i2, i3))
        end
        N = N2
    end

    # GeometryBasics.Mesh(GeometryBasics.meta(vertices; normals=ns), faces)
    return normal_mesh(bin_positions(B, radius), faces)
end


"""
    bin_positions(B::SSHBinner)

Returns a list of points corresponding to the center positions of each bin.
This list is equivalent to the vertices used in `vertex_mesh(B)`.

See also: [`vertex_mesh`](@ref)
"""
bin_positions(B::AbstractSSH, radius = 1.0) = bin_positions(B.tessellation, radius)
function bin_positions(B::SphereTessellationMap, radius)
    points = Vector{Point3f}(undef, length(B))
    points[1] = Point3f(0, 0, radius)

    idx = 2
    for i in 2:length(B.thetas)-2
        theta = 0.5(B.thetas[i] + B.thetas[i+1])
        N = B.phi_divisions[i]

        for k in 0:N-1
            phi = 2pi * (k + 0.5) / N
            points[idx] = to_cartesian(radius, theta, phi)
            idx += 1
        end
    end
    points[end] = Point3f(0, 0, -radius)

    return points
end


"""
    line_segments_minimal(B::AbstractSSH)

Generates linesegments marking each bin area. This version is minimal in the
sense that doesn't interpolate any points.

See also: [`line_segments`](@ref)
"""
line_segments_minimal(B::AbstractSSH, radius = 1) = line_segments_minimal(B.tessellation, radius)
function line_segments_minimal(B::SphereTessellationMap, radius)
    lines = Point3f[]

    # first ring
    theta = B.thetas[2]
    N = B.phi_divisions[2]
    for k in 0:N-1
        phi1 = 2pi * k / N
        phi2 = 2pi * (k+1) / N
        append!(lines, [to_cartesian(radius, theta, phi1), to_cartesian(radius, theta, phi2)])
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
            append!(lines, [to_cartesian(radius, theta2, phi1), to_cartesian(radius, theta2, phi2)])
            # vertical line
            append!(lines, [to_cartesian(radius, theta2, phi1), to_cartesian(radius, theta1, phi1)])
        end
    end

    return lines
end


"""
    line_segments(B::AbstractSSH; N_fragments = 32)

Generates line segments marking each bin area. Extra points may be interpolated
to generate smooth arcs.

Specifically, points are added to an arc, such that it would include at least
`N_fragments` points over an angle 0:pi.

See also: [`line_segments_minimal`](@ref)
"""
function line_segments(B::AbstractSSH, radius = 1.0; N_fragments=32)
    return line_segments(B.tessellation, radius, N_fragments)
end
function line_segments(B::SphereTessellationMap, radius, N_fragments)
    lines = Point3f[]

    # horizontal rings
    for theta in B.thetas[2:end-1]
        for k in 0:2N_fragments-1
            phi1 = pi * k / N_fragments
            phi2 = pi * (k+1) / N_fragments
            push!(lines, to_cartesian(radius, theta, phi1), to_cartesian(radius, theta, phi2))
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
            phi = 2pi * k / N
            for c in 0:substeps-1
                t1 = theta1 + c/substeps * (theta2 - theta1)
                t2 = theta1 + (c+1)/substeps * (theta2 - theta1)

                # vertical line
                push!(lines, to_cartesian(radius, t2, phi), to_cartesian(radius, t1, phi))
            end
        end
    end

    return lines
end
