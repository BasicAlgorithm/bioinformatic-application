module SECONDARY_STRUCTURE

using Dash
using DashHtmlComponents
using DashCoreComponents

using Colors

using Cairo, Fontconfig
using Compose
using Graphs, GraphPlot

using DashCytoscape

export create_matrix_sa
export score_and_directions_sa
export generate_secundary_structure

aa_colors = Dict([
  ("A", colorant"#c64d49"),
  ("G", colorant"#42d42e"),
  ("C", colorant"#e7da27"),
  ("T", colorant"#29dae0"),
  ("U", colorant"#2944e0"),
])

function create_matrix_sa(size_column)
  matrix_scores = Matrix{Int32}(undef, size_column, size_column)
  matrix_directions = Matrix{String}(undef, size_column, size_column)

  fill!(matrix_directions, "")
  fill!(matrix_scores, 0)

  return matrix_scores, matrix_directions
end # create_matrix

function getAlpha(a::Char, b::Char)
  if a == 'A' && b == 'U'
    return -1
  end
  if a == 'U' && b == 'A'
    return -1
  end

  if a == 'G' && b == 'C'
    return -1
  end
  if a == 'C' && b == 'G'
    return -1
  end

  return 0
end # getAlpha

function getV4(matrix_scores, f, c)
  minvalue = 1000
  # println("max v of int16 is ", minvalue)

  for x in f+1:c-1
    minvalue = min(matrix_scores[f, x] + matrix_scores[x+1, c], minvalue)
  end

  return minvalue
end # getV4

function score_and_directions_sa(matrix_scores, matrix_directions, nucleotides, sizes)
  f = 1
  c = 2
  while true
    t_f = f
    t_c = c

    while true
      # println(t_f,"-",t_c)
      # main operation
      v1 = matrix_scores[t_f, t_c-1]
      v2 = matrix_scores[t_f+1, t_c]
      v3 = matrix_scores[t_f+1, t_c-1] + getAlpha(nucleotides[1][t_c], nucleotides[1][t_f])
      v4 = getV4(matrix_scores, t_f, t_c)
      # println("=======")
      # println(t_f," ",t_c)
      # println("- ",nucleotides[1][t_c], " ", nucleotides[1][t_f]," alpha ",getAlpha(nucleotides[1][t_c], nucleotides[1][t_f]))
      # println(v1)
      # println(v2)
      # println(v3)
      # println(v4)

      val_min = min(v1, v2, v3, v4)

      # matrix directions
      if (val_min == v1)
        matrix_directions[t_f, t_c] *= "l"
      end
      if (val_min == v2)
        matrix_directions[t_f, t_c] *= "b"
      end
      if (val_min == v3)
        matrix_directions[t_f, t_c] *= "d"
      end
      if (val_min == v4)
        matrix_directions[t_f, t_c] *= "x"
      end

      #   # final assignment
      matrix_scores[t_f, t_c] = val_min


      t_f += 1
      t_c += 1

      if t_c > sizes
        break
      end
    end
    # println("---")
    c += 1

    if c > sizes
      break
    end
  end

  find_min = findmin(matrix_scores)
  value_min = find_min[1]
  value_min_pos = find_min[2]

  return findmin(matrix_scores)
end # score_and_directions

function recursiveAdd!(nucleotides, dynamic_plot, graph, matrix_directions, f, c)
  if matrix_directions[f, c] == ""
    return
  else
    if first(matrix_directions[f, c]) == 'l'
      recursiveAdd!(nucleotides, dynamic_plot, graph, matrix_directions, f, c - 1)
    elseif first(matrix_directions[f, c]) == 'b'
      recursiveAdd!(nucleotides, dynamic_plot, graph, matrix_directions, f + 1, c)
    elseif first(matrix_directions[f, c]) == 'd'
      add_edge!(graph, f, c)
      push!(dynamic_plot, Dict("data" => Dict("source" => string(nucleotides[1][f]) * string(f), "target" => string(nucleotides[1][c]) * string(c)), "classes" => "red"))
      recursiveAdd!(nucleotides, dynamic_plot, graph, matrix_directions, f + 1, c - 1)
    elseif first(matrix_directions[f, c]) == 'x'
      println("que hago si gana v4???")
    end
  end

end # recursiveAdd!

function generate_secundary_structure(matrix_directions, find_min, nucleotides)
  # dynamic plot
  dynamic_plot = []

  # println("haber ", nucleotides[1])
  graph = Graph(length(nucleotides[1]), 0)

  is_first = true

  nodo_labels = []
  for c in eachindex(nucleotides[1])
    push!(nodo_labels, string(nucleotides[1][c]))

    # println("=> id ", string(nucleotides[1][c]) * string(c))
    # println("=> fn ", string(nucleotides[1][c]))
    if c == 1
      push!(dynamic_plot, Dict("data" => Dict("id" => string(nucleotides[1][c]) * string(c), "firstname" => string(nucleotides[1][c]) * "."), "classes" => "text"))
    elseif c == length(nucleotides[1])
      push!(dynamic_plot, Dict("data" => Dict("id" => string(nucleotides[1][c]) * string(c), "firstname" => string(nucleotides[1][c]) * ".."), "classes" => "text"))
    else
      push!(dynamic_plot, Dict("data" => Dict("id" => string(nucleotides[1][c]) * string(c), "firstname" => string(nucleotides[1][c])), "classes" => "text"))
    end
  end

  # nodo_colors = []
  # for c in 1:nv(graph)
  #   push!(nodo_colors, get(aa_colors, nodo_labels[c], colorant"gray0"))
  # end

  for i in 1:nv(graph)-1
    add_edge!(graph, i, i + 1)
    # println("=> id ", string(nucleotides[1][c]) * string(c))
    # println("=> fn ", string(nucleotides[1][c]))
    push!(dynamic_plot, Dict("data" => Dict("source" => string(nucleotides[1][i]) * string(i), "target" => string(nucleotides[1][i+1]) * string(i + 1)), "classes" => "red"))
  end

  recursiveAdd!(nucleotides, dynamic_plot, graph, matrix_directions, find_min[2][1], find_min[2][2])

  # draw(PNG("./src/assets/plot.png", 16cm, 16cm),
  #   gplot(graph, nodelabel=nodo_labels, nodefillc=nodo_colors))

  # static plot
  # r = []
  # push!(r, html_div("Secondary Structure [Zhang H, et al. 2019]"))
  # push!(r, html_img(src="/assets/plot.png"))

  # dynamic plot
  dy_plot = []
  push!(dy_plot, cyto_cytoscape(
    layout=Dict("name" => "cose"),
    # grid random circle cose concentric preset breadthfirst
    style=Dict("width" => "100%", "height" => "400px"),
    elements=dynamic_plot,
    stylesheet=[
      Dict(
        "selector" => "node",
        "style" => Dict(
          "label" => "data(firstname)"
        )
      ),
      Dict(
        "selector" => ".text",
        "style" => Dict(
          "color" => "#FFFFFF",
        )
      ),
      Dict(
        "selector" => ".red",
        "style" => Dict(
          "color" => "red",
          "line-color" => "#FFFFFF"
        )
      ),
      Dict(
        "selector" => "[firstname *= 'A']",
        "style" => Dict(
          "background-color" => "#c64d49",
          # "shape" => "rectangle"
        )
      ),
      Dict(
        "selector" => "[firstname *= 'C']",
        "style" => Dict(
          "background-color" => "#e7da27",

          # "shape" => "rectangle"
        )
      ),
      Dict(
        "selector" => "[firstname *= 'G']",
        "style" => Dict(
          "background-color" => "#42d42e",
          # "shape" => "rectangle"
        )
      ),
      Dict(
        "selector" => "[firstname *= 'T']",
        "style" => Dict(
          "background-color" => "#29dae0",
          # "shape" => "rectangle"
        )
      ),
      Dict(
        "selector" => "[firstname *= 'U']",
        "style" => Dict(
          "background-color" => "#2944e0",
          # "shape" => "rectangle"
        )
      )
    ]
  ))

  return dy_plot
end # generate_secundary_structure

end