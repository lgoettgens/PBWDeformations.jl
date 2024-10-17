function symmetrize(f::FreeAssAlgElem)
    R = parent(f)
    g = zero(R)
    for (c, exp) in zip(AbstractAlgebra.coefficients(f), AbstractAlgebra.exponent_words(f))
        g += R(
            fill(c * QQ(1, factorial(length(exp))), factorial(length(exp))),
            [ind for ind in permutations(exp)],
        )

    end
    return g
end

function set_show_colorful_html(value::Bool)
    old_value = get_show_colorful_html()
    @set_preferences!("show_colorful_html" => value)
    return old_value
end

function get_show_colorful_html()
    return @load_preference("show_colorful_html", default = false)::Bool
end

Base.show(io::IO, ::MIME"text/html", A::ArcDiagram) = print(io, svg_string(A; colorful=get_show_colorful_html()))

function svg_path_string(v::ArcDiagramVertex, nv::ArcDiagramVertex, dims)
    (;margin, h, w) = dims
    x_v = (vertex_index(v) - 1) * w + margin
    x_nv = (vertex_index(nv) - 1) * w + margin
    y_v = is_upper_vertex(v) ? margin : h + margin
    y_nv = is_upper_vertex(nv) ? margin : h + margin

    function y_c(x1, x2)
        return (1-exp(-0.15*(1+abs(x1 - x2))))*0.85
    end

    if is_upper_vertex(v)
        if is_upper_vertex(nv)
            x_c1 = (8*x_v + 2*x_nv) / 10
            y_c1 = y_c(vertex_index(v), vertex_index(nv)) * h + margin
            x_c2 = (2*x_v + 8*x_nv) / 10
            y_c2 = y_c(vertex_index(v), vertex_index(nv)) * h + margin
            return "M $(x_v) $(y_v) C $(x_c1) $(y_c1), $(x_c2) $(y_c2), $(x_nv) $(y_nv)"
        else
            return "M $(x_v) $(y_v) L $(x_nv) $(y_nv)"
        end
    else
        if is_upper_vertex(nv)
            return "M $(x_v) $(y_v) L $(x_nv) $(y_nv)"
        else
            x_c1 = (8*x_v + 2*x_nv) / 10
            y_c1 = (1 - y_c(vertex_index(v), vertex_index(nv))) * h + margin
            x_c2 = (2*x_v + 8*x_nv) / 10
            y_c2 = (1 - y_c(vertex_index(v), vertex_index(nv))) * h + margin
            return "M $(x_v) $(y_v) C $(x_c1) $(y_c1), $(x_c2) $(y_c2), $(x_nv) $(y_nv)"
        end
    end
end

function svg_string(A::ArcDiagram; colorful=false)
    myblue = "#008"
    mygray = "#aaa"
    size = 40
    dims = (
        rad=0.1*size,
        margin=0.2*size,
        h=1*size,
        w=0.5*size,
        stroke_width=2, #px
    )

    height = dims.h + 2 * dims.margin
    width = (max(n_upper_vertices(A), n_lower_vertices(A)) - 1) * dims.w + 2 * dims.margin

    svg = "<svg height=\"$(height)\" width=\"$(width)\" style=\"vertical-align:-$(1/2 * height)px;\">"

    svg *= "<defs>"
    if A isa ArcDiagramDirected
        # marker to be used as an arrowhead
        svg *= """<marker id=\"arrow\" refX=\"2\" refY=\"2\" markerWidth=\"3\" markerHeight=\"4\" orient=\"auto\">
                    <path d=\"M 0 0 V 4 L 2 2 Z\" stroke=\"context-stroke\" fill=\"context-stroke\"/>
                  </marker>"""
    end
    svg *= "</defs>"

    i = 0
    for v in vertices(A)
        nv = neighbor(A, v)
        _vertex_lt(v, nv) || continue
        col = colorful ? "hsl($(105*i),100%,45%)" : myblue
        svg *= "<path d=\"$(svg_path_string(v, nv, dims))\" fill=\"transparent\" stroke=\"$col\" stroke-width=\"$(dims.stroke_width)px\" $(A isa ArcDiagramDirected ? "marker-end=\"url(#arrow)\"" : "")/>"
        i += 1
    end
    # for (i, n_verts) in enumerate((n_upper_vertices(A), n_lower_vertices(A)))
    #     for j in 1:n_verts
    #         col = colorful ? mygray : myblue
    #         svg *= "<circle cx=\"$(dims.margin+(j-1)*dims.w)\" cy=\"$(dims.margin+(i-1)*dims.h)\" r=\"$(dims.rad)\" fill=\"$col\" />"
    #     end
    # end
    svg *= "</svg>"
    return svg
end

function Base.show(io::IO, mime::MIME"text/html", x::Tuple{Int,ArcDiagram})
    print(io, "($(first(x)), ")
    show(io, mime, x[2])
    print(io, ")")
end
