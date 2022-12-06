module STAR_ALIGNMENT

using Dash
using DashHtmlComponents
using DashCoreComponents

export getMatricesScoreAlignments
export getDifferentPairsAlignments
export getFirstAlign
export getMSA

include("blosum.jl")
using .BLOSUM

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
  matrix_scores[1, 1] = get(biology_score, "**", 9)
  matrix_directions[1, 1] = "*"

  # set first row default values
  for a in 2:size_column

    matrix_scores[1, a] = matrix_scores[1, a-1] + get(biology_score, nucleotides[second][1] * nucleotides[first][a], 7)
    matrix_directions[1, a] *= "l"
  end

  # set first column default values
  for a in 2:size_row

    matrix_scores[a, 1] = matrix_scores[a-1, 1] + get(biology_score, nucleotides[second][a] * nucleotides[first][1], 7)
    matrix_directions[a, 1] *= "u"
  end

  return matrix_scores, matrix_directions
end # create_matrix

function score_and_directions!(matrix_scores, matrix_directions, nucleotides, sizes, i_c, i_f)
  for f = 2:sizes[i_f], c = 2:sizes[i_c]
    # left
    v1 = matrix_scores[f, c-1] - 2
    # diagonal
    v2 = if (nucleotides[i_c][c] == nucleotides[i_f][f])
      1
    else
      -1
    end
    v2 += matrix_scores[f-1, c-1]
    # up
    v3 = matrix_scores[f-1, c] - 2

    val_max = max(v1, v2, v3)

    # matrix directions
    if (val_max == v1)
      matrix_directions[f, c] *= "l"
    end
    if (val_max == v2)
      matrix_directions[f, c] *= "d"
    end
    if (val_max == v3)
      matrix_directions[f, c] *= "u"
    end

    # final assignment
    matrix_scores[f, c] = val_max
  end

  return matrix_scores[sizes[i_f], sizes[i_c]]
end # score_and_directions

function recursive_searching(f, c, seq1, seq2, matrix_directions, nucleotides, array_of_results, i_c, i_f)

  if f == 1 && c == 1
    push!(array_of_results, seq1)
    push!(array_of_results, seq2)
    return

  else
    if last(matrix_directions[f, c]) == 'u'
      recursive_searching(
        f - 1,
        c,
        "-" * seq1,
        nucleotides[i_f][f] * seq2,
        matrix_directions,
        nucleotides,
        array_of_results, i_c, i_f
      )
    end

    if last(matrix_directions[f, c]) == 'd'
      recursive_searching(
        f - 1,
        c - 1,
        nucleotides[i_c][c] * seq1,
        nucleotides[i_f][f] * seq2,
        matrix_directions,
        nucleotides,
        array_of_results, i_c, i_f
      )
    end

    if last(matrix_directions[f, c]) == 'l'
      recursive_searching(
        f,
        c - 1,
        nucleotides[i_c][c] * seq1,
        "-" * seq2,
        matrix_directions,
        nucleotides,
        array_of_results, i_c, i_f
      )
    end

  end # matrix_directions[f,c] == "*"
end # recursive_searching

function calculate_secuences_recursively(sizes, matrix_directions, nucleotides, c, f)
  array_of_results::Array{String} = []

  #using recursion
  recursive_searching(sizes[f], sizes[c], "", "", matrix_directions, nucleotides, array_of_results, c, f)

  # return value
  return array_of_results
end # calculate_secuences_recursively

function getMatricesScoreAlignments(nucleotides, sizes, kind_sequences)
  global_matrix_scores = Matrix{Integer}(undef, length(nucleotides), length(nucleotides))
  global_matrix_alignments = Matrix{Any}(undef, length(nucleotides), length(nucleotides))
  fill!(global_matrix_scores, 0)

  for c = 1:length(nucleotides), f = c+1:length(nucleotides)
    matrix_scores_t, matrix_directions_t = create_matrix(sizes[c], sizes[f], nucleotides, c, f, kind_sequences)
    score_t = score_and_directions!(matrix_scores_t, matrix_directions_t, nucleotides, sizes, c, f)
    global_matrix_scores[c, f] = score_t
    global_matrix_scores[f, c] = score_t
    array_of_results_t = calculate_secuences_recursively(sizes, matrix_directions_t, nucleotides, c, f)
    global_matrix_alignments[c, f] = array_of_results_t
    global_matrix_alignments[f, c] = [array_of_results_t[2], array_of_results_t[1]]
  end

  return global_matrix_scores, global_matrix_alignments

end

function getDifferentPairsAlignments(global_matrix_alignments, nucleotides, center)
  final_aligments = []
  star_alignments = []
  max_length_final_aligment = 0
  for index = 1:length(nucleotides)
    if index == center
      continue
    end
    if length(global_matrix_alignments[center, index][1]) > max_length_final_aligment
      max_length_final_aligment = length(global_matrix_alignments[center, index][1])
    end
    if length(global_matrix_alignments[center, index][2]) > max_length_final_aligment
      max_length_final_aligment = length(global_matrix_alignments[center, index][2])
    end
    push!(final_aligments, [global_matrix_alignments[center, index][1], global_matrix_alignments[center, index][2]])
    push!(star_alignments, html_table(children=[
      html_tr([html_td(i, className=i * "-color") for i in global_matrix_alignments[center, index][1]]),
      html_tr([html_td(i, className=i * "-color") for i in global_matrix_alignments[center, index][2]]),
    ]
    ))
  end
  return max_length_final_aligment, final_aligments, star_alignments
end

function getFirstAlign(final_aligments, max_length_final_aligment)
  first_align = ""
  for i in eachindex(final_aligments)
    # println(i, " | ", max_length_final_aligment, " - ", final_aligments[i][1])
    if max_length_final_aligment == length(final_aligments[i][1])
      if first_align == ""
        first_align = final_aligments[i][1]
      elseif final_aligments[i][1] != first_align
        i1 = 1
        i2 = 1
        tmp_align = ""
        while i1 < length(first_align) || i2 < length(final_aligments[i][1])
          if i1 > length(first_align)
            tmp_align *= "-"
            i2 += 1
          elseif i2 > length(final_aligments[i][1])
            tmp_align *= "-"
            i1 += 1
          elseif first_align[i1] == final_aligments[i][1][i2]
            tmp_align *= first_align[i1]
            i1 += 1
            i2 += 1
          elseif first_align[i1] == '-'
            tmp_align *= "-"
            i1 += 1
          elseif final_aligments[i][1][i2] == '-'
            tmp_align *= "-"
            i2 += 1
          end
        end
        # for j in 1:length(first_align)+1
        #   if first_align[i1] == final_aligments[i][1][i2]
        #     tmp_align *= first_align[i1]
        #     i1 += 1
        #     i2 += 1
        #   elseif first_align[i1] == '-'
        #     tmp_align *= "-"
        #     i1 += 1
        #   elseif final_aligments[i][1][i2] == '-'
        #     tmp_align *= "-"
        #     i2 += 1
        #   end
        # end
        first_align = tmp_align
      end
    end
  end
  return first_align
end

function getMSA(final_aligments, first_align)
  msa = []
  # push!(msa, html_div(first_align))
  push!(msa, html_tr([html_td(i, className=i * "-color") for i in first_align]))

  for i in eachindex(final_aligments)
    tmp_align2 = ""
    it = 1
    for j in first_align

      if it > length(final_aligments[i][1])
        tmp_align2 *= "-"
        continue
      end

      if j == final_aligments[i][1][it]
        tmp_align2 *= string(final_aligments[i][2][it])
        it += 1
      else
        tmp_align2 *= "-"
      end

    end
    # push!(msa, html_div(tmp_align2))
    push!(msa, html_tr([html_td(i, className=i * "-color") for i in tmp_align2]))
  end
  return msa
end

end