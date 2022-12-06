module GLOBAL_LOCAL_ALIGNMENT

using Dash
using DashHtmlComponents
using DashCoreComponents

include("blosum.jl")
using .BLOSUM

export open_save_data_file1
export open_save_data_file
export create_matrix
export score_and_directions
export recursive_searching2
export calculate_secuences_recurvibily
export getIdentity
export multiColorAlignments
export create_matrix_la
export scoring_local_alignment
export recursive_searching
export calculate_matches

const global match::Int16 = 1
const global mismatch::Int16 = -1
const global diff::Int16 = -2
global q_results = 0

function open_save_data_file1(chain01, chain02)
  sizes = [length(chain01) + 1, length(chain02) + 1]
  nucleotides = ["*" * chain01, "*" * chain02]

  # println(sizes)
  # println(nucleotides)

  return sizes, nucleotides
end # open_save_data_file

function open_save_data_file(chain01, chain02, chain03)
  sizes = [length(chain01) + 1, length(chain02) + 1, length(chain03) + 1]
  nucleotides = ["*" * chain01, "*" * chain02, "*" * chain03]

  # println(sizes)
  # println(nucleotides)

  return sizes, nucleotides
end # open_save_data_file

function create_matrix(size_column, size_row, nucleotides, first, second, kind_sequences)
  biology_score = nothing
  if kind_sequences == "PROTEIN"
    biology_score = blosum62
  else
    biology_score = defaultAdn
  end

  # println("\n DEBUG SCORE SELECTED ", get(biology_score, "AA", 9))

  matrix_scores = Matrix{Float32}(undef, size_row, size_column)
  matrix_directions = Matrix{String}(undef, size_row, size_column)
  fill!(matrix_directions, "")

  #set predefined values to start calculating
  matrix_scores[1, 1] = get(biology_score, "**", 99)
  matrix_directions[1, 1] = "*"

  # set first row default values
  for a in 2:size_column

    matrix_scores[1, a] = matrix_scores[1, a-1] + get(biology_score, nucleotides[second][1] * nucleotides[first][a], 77)
    matrix_directions[1, a] *= "l"
  end

  # set first column default values
  for a in 2:size_row

    matrix_scores[a, 1] = matrix_scores[a-1, 1] + get(biology_score, nucleotides[second][a] * nucleotides[first][1], 77)
    matrix_directions[a, 1] *= "u"
  end

  return matrix_scores, matrix_directions
end # create_matrix

function score_and_directions(matrix_scores, matrix_directions, nucleotides, sizes, kind_sequences)
  biology_score = nothing
  if kind_sequences == "PROTEIN"
    biology_score = blosum62
  else
    biology_score = defaultAdn
  end

  for f = 2:sizes[2], c = 2:sizes[1]
    # left
    v1 = matrix_scores[f, c-1] + get(biology_score, nucleotides[2][f] * "*", 99)
    # diagonal
    v2 = matrix_scores[f-1, c-1] + get(biology_score, nucleotides[2][f] * nucleotides[1][c], 99)
    # up
    v3 = matrix_scores[f-1, c] + get(biology_score, nucleotides[1][c] * "*", 99)

    val_min = max(v1, v2, v3)

    # matrix directions
    if (val_min == v1)
      matrix_directions[f, c] *= "l"
    end
    if (val_min == v2)
      matrix_directions[f, c] *= "d"
    end
    if (val_min == v3)
      matrix_directions[f, c] *= "u"
    end

    # final assignment
    matrix_scores[f, c] = val_min
  end
end # score_and_directions

function recursive_searching2(f, c, seq1, seq2, matrix_directions, nucleotides, array_of_results, max_results)

  if f == 1 && c == 1
    push!(array_of_results, seq1)
    push!(array_of_results, seq2)
    return

  else

    for nucleotide in matrix_directions[f, c]
      if nucleotide == 'u'
        recursive_searching2(
          f - 1,
          c,
          "-" * seq1,
          nucleotides[2][f] * seq2,
          matrix_directions,
          nucleotides,
          array_of_results, max_results
        )
      end

      if length(array_of_results) >= max_results * 2
        return
      end

      if nucleotide == 'd'
        recursive_searching2(
          f - 1,
          c - 1,
          nucleotides[1][c] * seq1,
          nucleotides[2][f] * seq2,
          matrix_directions,
          nucleotides,
          array_of_results, max_results
        )
      end

      if length(array_of_results) >= max_results * 2
        return
      end

      if nucleotide == 'l'
        recursive_searching2(
          f,
          c - 1,
          nucleotides[1][c] * seq1,
          "-" * seq2,
          matrix_directions,
          nucleotides,
          array_of_results, max_results
        )
      end

      if length(array_of_results) >= max_results * 2
        return
      end

    end
  end # matrix_directions[f,c] == "*"
end # recursive_searching

function calculate_secuences_recurvibily(sizes, matrix_directions, nucleotides, max_results)
  array_of_results::Array{String} = []

  #using recursion
  recursive_searching2(sizes[2], sizes[1], "", "", matrix_directions, nucleotides, array_of_results, max_results)

  # return value
  return array_of_results
end # calculate_secuences_recurvibily

function getIdentity(a)
  if a == ""
    # println("sequence empty")
    return ""
  end
  # protein = false
  adn = false
  arn = false
  # println(a)
  # println(typeof(a))
  for c in a
    c = string(c)
    # println(c)
    # println(typeof(c))
    if findnext(c, "RNDQEHILKMFPSWYVBZX", 1) !== nothing
      return "PROTEIN"
    elseif findnext(c, "U", 1) !== nothing
      arn = true
    end
  end

  if arn == true
    return "RNA"
  else
    return "DNA"
  end

end # identified

function multiColorAlignments(array_of_results)
  children = []
  first_color = true
  q_results = trunc(Int, length(array_of_results) / 2 - 1)
  # push!(children, html_div("Quantity is: $(q_results + 1)"))
  push!(children, html_br())

  for c = 0:q_results
    kids = []
    for y in array_of_results[c*2+1]
      if first_color
        push!(kids, html_span(y, className="F-color"))
        first_color = false
      else
        push!(kids, html_span(y, className="S-color"))
        first_color = true
      end
    end
    push!(kids, html_br())

    kids2 = []
    first_color = true
    for y in array_of_results[c*2+2]
      if first_color
        push!(kids2, html_span(y, className="F-color"))
        first_color = false
      else
        push!(kids2, html_span(y, className="S-color"))
        first_color = true
      end
    end
    push!(kids2, html_br())
    push!(kids2, html_br())

    kidsf = []
    append!(kidsf, kids)
    append!(kidsf, kids2)

    # push!(children, html_div(children=kidsf))
    push!(children, html_table(children=[
      html_tr([html_td(i) for i in kids]),
      html_tr([html_td(i) for i in kids2]),
    ]))
  end

  return children
end # multiColorAlignments

# LOCAL ALIGNMENT 
function create_matrix_la(size_column, size_row)
  # matrix_scores = fill(size_column * size_row * -1, (size_row, size_column))
  matrix_scores = Matrix{Int32}(undef, size_row, size_column)
  # delete matrix for free ram memory
  # matrix_scores = nothing
  # GC.gc()
  matrix_directions = Matrix{String}(undef, size_row, size_column)
  fill!(matrix_directions, "")
  fill!(matrix_scores, 0)

  #set predefined values to start calculating
  matrix_directions[1, 1] = "*"

  return matrix_scores, matrix_directions
end # create_matrix

function scoring_local_alignment(matrix_scores, matrix_directions, nucleotides, sizes, kind_sequences)
  biology_score = defaultAdn

  # if kind_sequences == "PROTEIN"
  #   biology_score = blosum62
  # else
  #   biology_score = defaultAdn
  # end

  max_value = 0
  max_value_array = []

  for f = 2:sizes[2], c = 2:sizes[1]
    # left
    v1 = matrix_scores[f, c-1] + get(biology_score, nucleotides[2][f] * "*", -2)
    # diagonal
    # v2 = matrix_scores[f-1, c-1] + get(biology_score, nucleotides[2][f] * nucleotides[1][c], 99)
    v2 = if (nucleotides[1][c] == nucleotides[2][f])
      1
    else
      -1
    end
    v2 += matrix_scores[f-1, c-1]
    # up
    v3 = matrix_scores[f-1, c] + get(biology_score, nucleotides[1][c] * "*", -2)

    val_min = max(0, v1, v2, v3)

    # matrix directions
    if (val_min == v1)
      matrix_directions[f, c] *= "l"
    end
    if (val_min == v2)
      matrix_directions[f, c] *= "d"
    end
    if (val_min == v3)
      matrix_directions[f, c] *= "u"
    end

    # final assignment
    matrix_scores[f, c] = val_min

    # fill returns
    if val_min > max_value
      max_value_array = []
      max_value = val_min
      push!(max_value_array, [f, c])
    elseif val_min == max_value
      push!(max_value_array, [f, c])
    else

    end

  end
  return max_value, max_value_array
end # scoring

function recursive_searching(f, limit, nucleotides)
  if f == limit
    return ""
  else
    return recursive_searching(f - 1, limit, nucleotides) * nucleotides[2][f]
  end
end # recursive_searching

function calculate_matches(max_value, max_value_array, nucleotides)
  array_of_results::Array{String} = []
  r = []

  for pair_index in max_value_array
    push!(array_of_results, recursive_searching(pair_index[1], pair_index[1] - max_value, nucleotides))
  end

  push!(r, html_div("There are " * string(length(array_of_results)) * " start points with max value " * string(max_value)))
  push!(r, html_br())
  for i in array_of_results
    push!(r, html_div(i, className="content"))
  end

  return r
end # calculate_matches

end