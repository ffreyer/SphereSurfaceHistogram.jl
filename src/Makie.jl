@reexport using GeometryTypes, Colors, Makie


function to_rects(B::SSHBinner, extrude=0)
    thetas = B.thetas
    phi_divs = B.phi_divisions
    dphis = 2pi ./ phi_divs

    vertices = Point3f0[]
    faces = Face{3, Int}[]

    push!(vertices, Point3f0(0, 0, 1) + extrude * Point3f0(0,0,1))
    theta2 = thetas[2]
    for k in 0:phi_divs[2]-1
        phi1 = dphis[2] * k
        phi2 = dphis[2] * (k+1)
        pos_1 = to_cartesian(theta2, phi1) + extrude * Point3f0(0,0,1)
        pos_2 = to_cartesian(theta2, phi2) + extrude * Point3f0(0,0,1)
        push!(vertices, pos_1)
        push!(vertices, pos_2)
        push!(faces, Face(1, 2+2k, 3+2k))
    end

    for i in 2:length(dphis)-1
        theta1 = thetas[i]
        theta2 = thetas[i+1]
        for k in 0:phi_divs[i]-1
            phi1 = dphis[i] * k
            phi2 = dphis[i] * (k+1)
            pos_11 = to_cartesian(theta1, phi1)
            pos_12 = to_cartesian(theta1, phi2)
            pos_21 = to_cartesian(theta2, phi1)
            pos_22 = to_cartesian(theta2, phi2)
            # delete me if you want
            dir = 0.25(pos_11 + pos_12 + pos_21 + pos_22)
            pos_11 += extrude*dir
            pos_12 += extrude*dir
            pos_21 += extrude*dir
            pos_22 += extrude*dir
            # /delete me if you want
            l = length(vertices)
            append!(vertices, [pos_11, pos_12, pos_21, pos_22])
            push!(faces, Face(l+1, l+3, l+2))
            push!(faces, Face(l+2, l+3, l+4))
        end
    end

    push!(vertices, Point3f0(0, 0, -1) + extrude * Point3f0(0,0,-1))
    l = length(vertices)
    theta2 = thetas[end-1]
    for k in 0:phi_divs[end-1]-1
        phi1 = dphis[2] * k
        phi2 = dphis[2] * (k+1)
        pos_1 = to_cartesian(theta2, phi1) + extrude * Point3f0(0,0,-1)
        pos_2 = to_cartesian(theta2, phi2) + extrude * Point3f0(0,0,-1)
        push!(vertices, pos_1)
        push!(vertices, pos_2)
        push!(faces, Face(l, l+2+2k, l+1+2k))
    end

    GLNormalMesh(vertices = vertices, faces = faces)
end


function to_dual_mesh(B::SSHBinner)
    thetas = B.thetas
    phi_divs = B.phi_divisions
    dphis = 2pi ./ phi_divs

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
    theta = 0.5(thetas[2] + thetas[3])
    for k in 0:phi_divs[2]-1
        phi = dphis[2] * (k + 0.5)
        push!(vertices, to_cartesian(theta, phi))
        push!(faces, Face(1, k+2, ((k+1)%phi_divs[2])+2))
    end

    for i in 3:length(thetas)-2
        theta = 0.5(thetas[i] + thetas[i+1])

        for k in 0:phi_divs[i]-1
            phi = dphis[i] * (k + 0.5)
            push!(vertices, to_cartesian(theta, phi))

            l = length(vertices)
            l < 3 && continue
            if phi_divs[i-1] == phi_divs[i]
                push!(faces, Face(l-1, l, l-phi_divs[i]))
                push!(faces, Face(l, l-phi_divs[i]+1, l-phi_divs[i]))

            elseif phi_divs[i-1] < phi_divs[i]
                # offset to top layer
                # if factor = 2: increment inner_index every 2nd step
                # if factor = 4: increment inner_index every 4th step
                factor = div(phi_divs[i], phi_divs[i-1])
                # equivalent to l - phi_divs - div(factor-1 + k*(factor-1), factor)
                inner_index = l - phi_divs[i-1] - div((k+1)*(factor-1), factor)
                # println("l = $l   \t  inner_index = $inner_index")
                if phi_divs[i]-1 > k > 0
                    if k % 2 == 1
                        push!(faces, Face(l-1, l, inner_index))
                        # push!(faces, Face(l, l-1, inner_index))
                        # push!(faces, Face(l, inner_index, inner_index+1))
                        push!(faces, Face(l, inner_index+1, inner_index))
                    else
                        push!(faces, Face(l-1, l, inner_index))
                    end
                elseif k == 0
                    last_index = l + phi_divs[i] - 1
                    push!(faces, Face(last_index-1, last_index, l-1))
                    push!(faces, Face(last_index, inner_index, l-1))
                    push!(faces, Face(last_index, l, inner_index))
                end

            elseif phi_divs[i-1] > phi_divs[i]
                offset = phi_divs[i-1] - k
                if k > 0
                    push!(faces, Face(l-1, l, l-offset-1))
                    push!(faces, Face(l, l-offset, l-offset-1))
                    push!(faces, Face(l, l-offset+1, l-offset))
                else
                    push!(faces, Face(l-1, l, l-offset))
                    push!(faces, Face(l, l-offset, l-offset+1))
                    push!(faces, Face(l, l-1, l+phi_divs[i]-1))
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
    for k in 0:phi_divs[end-1]-1
        push!(faces, Face(l, l-k-1, l - ((k+1)%phi_divs[end-1]) - 1))
    end

    GLNormalMesh(vertices = vertices, faces = faces)
end


function dual_points(B::SSHBinner)
    thetas = B.thetas
    phi_divs = B.phi_divisions
    dphis = 2pi ./ phi_divs

    points = Point3f0[Point3f0(0, 0, 1)]

    for i in 2:length(thetas)-2
        theta = 0.5(thetas[i] + thetas[i+1])
        for k in 0:phi_divs[i]-1
            phi = dphis[i] * (k + 0.5)
            push!(points, to_cartesian(theta, phi))
        end
    end
    push!(points, Point3f0(0, 0, -1))

    points
end
