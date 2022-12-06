using Dash
using DashHtmlComponents
using DashCoreComponents

using Random

include("global_local_alignment.jl")
using .GLOBAL_LOCAL_ALIGNMENT

include("star_alignment.jl")
using .STAR_ALIGNMENT

include("secondary_structure.jl")
using .SECONDARY_STRUCTURE

include("utilities.jl")
using .UTILITIES

########## WEBAPP setup ##########
app = dash(suppress_callback_exceptions=true)
app.title = "Bioinformatic App"

########## WEBAPP first layout ##########
app.layout = html_div(style=Dict("width" => "100%")) do
    # html_img(src="/assets/logo.png"),
    html_div("Bioinformatic Application", className="super-title"),
    html_div([
            dcc_dropdown(id="choosen_feature", options=dropdown_initial_options, placeholder="choose type of analysis", value="")
        ], className="dropdown"),
    html_div(id="features_app"),
    html_footer("Miguel Nieto | 2022", className="footer")
end

########## WEBAPP callbacks ##########

# general callback
callback!(
    app,
    Output("features_app", "children"),
    Input("choosen_feature", "value"),
    prevent_initial_call=true,
) do option_selected
    return get(functions_tags, string(option_selected), html_div("error"))
end

# Callback for Global and local alignments
callback!(
    app,
    Output("lenght_01", "children"),
    Output("lenght_02", "children"),
    Output("identified_01", "children"),
    Output("identified_02", "children"),
    Output("alignments", "children"),
    Output("score", "children"),
    Output("la_aligments", "children"),
    Input("button_global_local_alig", "n_clicks"),
    State("input_01", "value"),
    State("input_02", "value"),
    State("slider-value", "value"),
    prevent_initial_call=true,) do n_clicks, chain01, chain02, max_alignments
    id_chain_1 = getIdentity(uppercase(chain01))
    id_chain_2 = getIdentity(uppercase(chain02))
    kind_sequences = ""
    if id_chain_1 == "PROTEIN" || id_chain_2 == "PROTEIN"
        kind_sequences = "PROTEIN"
    elseif id_chain_1 == "RNA" || id_chain_2 == "RNA"
        kind_sequences = "RNA"
    else
        kind_sequences = "DNA"
    end
    ################## GLOBAL ALIGNMENT
    sizes, nucleotides = open_save_data_file1(uppercase(chain01), uppercase(chain02))
    matrix_scores, matrix_directions = create_matrix(sizes[1], sizes[2], nucleotides, 1, 2, kind_sequences)
    score_and_directions(matrix_scores, matrix_directions, nucleotides, sizes, kind_sequences)
    array_of_results = calculate_secuences_recurvibily(sizes, matrix_directions, nucleotides, max_alignments)
    alignments = multiColorAlignments(array_of_results)

    # show(stdout, "text/plain", matrix_scores)
    # show(stdout, "text/plain", matrix_directions)

    ################## LOCAL ALIGNMENT
    matrix_scores_la, matrix_directions_la = create_matrix_la(sizes[1], sizes[2])
    max_value, max_value_array = scoring_local_alignment(matrix_scores_la, matrix_directions_la, nucleotides, sizes, kind_sequences)

    # show(stdout, "text/plain", matrix_scores_la)
    # show(stdout, "text/plain", matrix_directions_la)
    # show(stdout, "text/plain", max_value)
    # println("\nDEBUG max value ", max_value)
    # show(stdout, "text/plain", max_value_array)

    array_of_results_la = calculate_matches(max_value, max_value_array, nucleotides)

    return (length(chain01), length(chain02),
        id_chain_1, id_chain_2,
        alignments, matrix_scores[sizes[2], sizes[1]], array_of_results_la,
    )
end # callback!

# Callback for create multi input sequence
callback!(
    app,
    Output("msa_with_multi_sequence", "children"),
    Input("q_sequence_msa", "value"),
    prevent_initial_call=true,
) do quantity_selected
    inp = []
    child_identifiers = []
    i = 1
    for n in quantity_selected*10:quantity_selected*10+quantity_selected-1
        push!(inp, html_div(
            children=[
                dcc_input(id="input_" * string(n), type="text", value="ACCCGTU", className="form__field"),
                html_label("Write your " * string(i) * " sequence", className="form__label"),
            ], className="form__group field")
        )
        i += 1
    end

    push!(inp, html_div(
        children=[
            html_button(id="button_star_alignment"* string(quantity_selected), children="Align sequences", className="button-89", n_clicks=0),
        ], className="div_button"))
    push!(inp, html_br())
    push!(inp, html_br())
    push!(inp, html_br())
    push!(inp, html_div(id="card_output_" * string(quantity_selected), className="mini-cards"))
    push!(inp, html_div(
        html_h1("Alignment choosen"), className="card-title"))
    push!(inp, html_div(
        children=[
            html_div(id="sa_align_choosen"* string(quantity_selected)),
        ], className="card"))
    push!(inp, html_div(
        html_h1("Pair alignments individually"), className="card-title"))
    push!(inp, html_div(
        children=[
            html_div(id="sa_pairs_aligments"* string(quantity_selected)),
        ], className="card"))
    push!(inp, html_div(
        html_h1("Multiple alignment MSA"), className="card-title"))
    push!(inp, html_div(
        children=[
            html_div(id="sa_multiple_aligments"* string(quantity_selected)),
        ], className="card"))

    return inp
end

# Callback for 5 Star alignments
callback!(
    app,
    Output("card_output_5", "children"),
    Output("sa_align_choosen5", "children"),
    Output("sa_pairs_aligments5", "children"),
    Output("sa_multiple_aligments5", "children"),
    Input("button_star_alignment5", "n_clicks"),
    State("input_50", "value"),
    State("input_51", "value"),
    State("input_52", "value"),
    State("input_53", "value"),
    State("input_54", "value"),
    prevent_initial_call=true,) do n_clicks, chain01, chain02, chain03, chain04, chain05
    sizes = []
    nucleotides = []
    output_cards = []
    if chain01 != ""
        # println("es string vacio")
        push!(nucleotides, "*" * uppercase(chain01))
        push!(sizes, length(chain01) + 1)
        push!(output_cards, html_div(
            children=[
                html_div("First sequence type", className="title"),
                html_br(),
                html_div(getIdentity(uppercase(chain01)), className="content"),
            ], className="mini-card")
        )
    end
    if chain02 != ""
        # println("es string vacio")
        push!(nucleotides, "*" * uppercase(chain02))
        push!(sizes, length(chain02) + 1)
        push!(output_cards, html_div(
            children=[
                html_div("Second sequence type", className="title"),
                html_br(),
                html_div(getIdentity(uppercase(chain02)), className="content"),
            ], className="mini-card")
        )

    end
    if chain03 != ""
        # println("es string vacio")
        push!(nucleotides, "*" * uppercase(chain03))
        push!(sizes, length(chain03) + 1)
        push!(output_cards, html_div(
            children=[
                html_div("Third sequence type", className="title"),
                html_br(),
                html_div(getIdentity(uppercase(chain03)), className="content"),
            ], className="mini-card")
        )

    end
    if chain04 != ""
        # println("es string vacio")
        push!(nucleotides, "*" * uppercase(chain04))
        push!(sizes, length(chain04) + 1)
        push!(output_cards, html_div(
            children=[
                html_div("Fourth sequence type", className="title"),
                html_br(),
                html_div(getIdentity(uppercase(chain04)), className="content"),
            ], className="mini-card")
        )

    end
    if chain05 != ""
        # println("es string vacio")
        push!(nucleotides, "*" * uppercase(chain05))
        push!(sizes, length(chain05) + 1)
        push!(output_cards, html_div(
            children=[
                html_div("Fifth sequence type", className="title"),
                html_br(),
                html_div(getIdentity(uppercase(chain05)), className="content"),
            ], className="mini-card")
        )

    end

    kind_sequences = ""
    if getIdentity(uppercase(chain01)) == "PROTEIN" || getIdentity(uppercase(chain02)) == "PROTEIN" || getIdentity(uppercase(chain03)) == "PROTEIN" || getIdentity(uppercase(chain04)) == "PROTEIN" || getIdentity(uppercase(chain05)) == "PROTEIN"
        kind_sequences = "PROTEIN"
    elseif getIdentity(uppercase(chain01)) == "RNA" || getIdentity(uppercase(chain02)) == "RNA" || getIdentity(uppercase(chain03)) == "RNA" || getIdentity(uppercase(chain04)) == "RNA" || getIdentity(uppercase(chain05)) == "RNA"
        kind_sequences = "RNA"
    else
        kind_sequences = "DNA"
    end

    push!(output_cards, html_div(
            children=[
                html_div("SEQUENCES ALIGNED AS", className="title"),
                html_br(),
                html_div(kind_sequences, className="content"),
            ], className="mini-card")
        )

    global_matrix_scores, global_matrix_alignments = getMatricesScoreAlignments(nucleotides, sizes, kind_sequences)

    mm = findmax(sum(global_matrix_scores, dims=2))
    alignment_choose = string(mm[2][1]) * "th sequence choosen => " * nucleotides[mm[2][1]][begin+1:end] * " => with score " * string(mm[1])

    max_length_final_aligment, final_aligments, pairs_alignments_tags = getDifferentPairsAlignments(global_matrix_alignments, nucleotides, mm[2][1])

    # creating first and base alignment
    first_align = getFirstAlign(final_aligments, max_length_final_aligment)

    # creating msa
    msa = getMSA(final_aligments, first_align)

    # show(stdout, "text/plain", msa)
    return (output_cards, alignment_choose, pairs_alignments_tags, html_table(children=msa))
end # callback!

# Callback for 4 Star alignments
callback!(
    app,
    Output("card_output_4", "children"),
    Output("sa_align_choosen4", "children"),
    Output("sa_pairs_aligments4", "children"),
    Output("sa_multiple_aligments4", "children"),
    Input("button_star_alignment4", "n_clicks"),
    State("input_40", "value"),
    State("input_41", "value"),
    State("input_42", "value"),
    State("input_43", "value"),
    prevent_initial_call=true,) do n_clicks, chain01, chain02, chain03, chain04
    sizes = []
    nucleotides = []
    output_cards = []
    if chain01 != ""
        # println("es string vacio")
        push!(nucleotides, "*" * uppercase(chain01))
        push!(sizes, length(chain01) + 1)
        push!(output_cards, html_div(
            children=[
                html_div("First sequence type", className="title"),
                html_br(),
                html_div(getIdentity(uppercase(chain01)), className="content"),
            ], className="mini-card")
        )
    end
    if chain02 != ""
        # println("es string vacio")
        push!(nucleotides, "*" * uppercase(chain02))
        push!(sizes, length(chain02) + 1)
        push!(output_cards, html_div(
            children=[
                html_div("Second sequence type", className="title"),
                html_br(),
                html_div(getIdentity(uppercase(chain02)), className="content"),
            ], className="mini-card")
        )

    end
    if chain03 != ""
        # println("es string vacio")
        push!(nucleotides, "*" * uppercase(chain03))
        push!(sizes, length(chain03) + 1)
        push!(output_cards, html_div(
            children=[
                html_div("Third sequence type", className="title"),
                html_br(),
                html_div(getIdentity(uppercase(chain03)), className="content"),
            ], className="mini-card")
        )

    end
    if chain04 != ""
        # println("es string vacio")
        push!(nucleotides, "*" * uppercase(chain04))
        push!(sizes, length(chain04) + 1)
        push!(output_cards, html_div(
            children=[
                html_div("Fourth sequence type", className="title"),
                html_br(),
                html_div(getIdentity(uppercase(chain04)), className="content"),
            ], className="mini-card")
        )

    end

    kind_sequences = ""
    if getIdentity(uppercase(chain01)) == "PROTEIN" || getIdentity(uppercase(chain02)) == "PROTEIN" || getIdentity(uppercase(chain03)) == "PROTEIN" || getIdentity(uppercase(chain04)) == "PROTEIN" 
        kind_sequences = "PROTEIN"
    elseif getIdentity(uppercase(chain01)) == "RNA" || getIdentity(uppercase(chain02)) == "RNA" || getIdentity(uppercase(chain03)) == "RNA" || getIdentity(uppercase(chain04)) == "RNA"
        kind_sequences = "RNA"
    else
        kind_sequences = "DNA"
    end

    push!(output_cards, html_div(
            children=[
                html_div("SEQUENCES ALIGNED AS", className="title"),
                html_br(),
                html_div(kind_sequences, className="content"),
            ], className="mini-card")
        )

    global_matrix_scores, global_matrix_alignments = getMatricesScoreAlignments(nucleotides, sizes, kind_sequences)

    mm = findmax(sum(global_matrix_scores, dims=2))
    alignment_choose = string(mm[2][1]) * "th sequence choosen => " * nucleotides[mm[2][1]][begin+1:end] * " => with score " * string(mm[1])

    max_length_final_aligment, final_aligments, pairs_alignments_tags = getDifferentPairsAlignments(global_matrix_alignments, nucleotides, mm[2][1])

    # creating first and base alignment
    first_align = getFirstAlign(final_aligments, max_length_final_aligment)

    # creating msa
    msa = getMSA(final_aligments, first_align)

    # show(stdout, "text/plain", msa)
    return (output_cards, alignment_choose, pairs_alignments_tags, html_table(children=msa))
end # callback!

# Callback for 3 Star alignments
callback!(
    app,
    Output("card_output_3", "children"),
    Output("sa_align_choosen3", "children"),
    Output("sa_pairs_aligments3", "children"),
    Output("sa_multiple_aligments3", "children"),
    Input("button_star_alignment3", "n_clicks"),
    State("input_30", "value"),
    State("input_31", "value"),
    State("input_32", "value"),
    prevent_initial_call=true,) do n_clicks, chain01, chain02, chain03
    sizes = []
    nucleotides = []
    output_cards = []
    if chain01 != ""
        # println("es string vacio")
        push!(nucleotides, "*" * uppercase(chain01))
        push!(sizes, length(chain01) + 1)
        push!(output_cards, html_div(
            children=[
                html_div("First sequence type", className="title"),
                html_br(),
                html_div(getIdentity(uppercase(chain01)), className="content"),
            ], className="mini-card")
        )
    end
    if chain02 != ""
        # println("es string vacio")
        push!(nucleotides, "*" * uppercase(chain02))
        push!(sizes, length(chain02) + 1)
        push!(output_cards, html_div(
            children=[
                html_div("Second sequence type", className="title"),
                html_br(),
                html_div(getIdentity(uppercase(chain02)), className="content"),
            ], className="mini-card")
        )
    end
    if chain03 != ""
        # println("es string vacio")
        push!(nucleotides, "*" * uppercase(chain03))
        push!(sizes, length(chain03) + 1)
        push!(output_cards, html_div(
            children=[
                html_div("Third sequence type", className="title"),
                html_br(),
                html_div(getIdentity(uppercase(chain03)), className="content"),
            ], className="mini-card")
        )
    end


    kind_sequences = ""
    if getIdentity(uppercase(chain01)) == "PROTEIN" || getIdentity(uppercase(chain02)) == "PROTEIN" || getIdentity(uppercase(chain03)) == "PROTEIN" 
        kind_sequences = "PROTEIN"
    elseif getIdentity(uppercase(chain01)) == "RNA" || getIdentity(uppercase(chain02)) == "RNA" || getIdentity(uppercase(chain03)) == "RNA"
        kind_sequences = "RNA"
    else
        kind_sequences = "DNA"
    end

    push!(output_cards, html_div(
            children=[
                html_div("SEQUENCES ALIGNED AS", className="title"),
                html_br(),
                html_div(kind_sequences, className="content"),
            ], className="mini-card")
        )

    global_matrix_scores, global_matrix_alignments = getMatricesScoreAlignments(nucleotides, sizes, kind_sequences)

    mm = findmax(sum(global_matrix_scores, dims=2))
    alignment_choose = string(mm[2][1]) * "th sequence choosen => " * nucleotides[mm[2][1]][begin+1:end] * " => with score " * string(mm[1])

    max_length_final_aligment, final_aligments, pairs_alignments_tags = getDifferentPairsAlignments(global_matrix_alignments, nucleotides, mm[2][1])

    # creating first and base alignment
    first_align = getFirstAlign(final_aligments, max_length_final_aligment)

    # creating msa
    msa = getMSA(final_aligments, first_align)

    # show(stdout, "text/plain", msa)
    return (output_cards, alignment_choose, pairs_alignments_tags, html_table(children=msa))
end # callback!

# Callback for 2 Star alignments
callback!(
    app,
    Output("card_output_2", "children"),
    Output("sa_align_choosen2", "children"),
    Output("sa_pairs_aligments2", "children"),
    Output("sa_multiple_aligments2", "children"),
    Input("button_star_alignment2", "n_clicks"),
    State("input_20", "value"),
    State("input_21", "value"),
    prevent_initial_call=true,) do n_clicks, chain01, chain02
    sizes = []
    nucleotides = []
    output_cards = []
    if chain01 != ""
        # println("es string vacio")
        push!(nucleotides, "*" * uppercase(chain01))
        push!(sizes, length(chain01) + 1)
        push!(output_cards, html_div(
            children=[
                html_div("First sequence type", className="title"),
                html_br(),
                html_div(getIdentity(uppercase(chain01)), className="content"),
            ], className="mini-card")
        )
    end
    if chain02 != ""
        # println("es string vacio")
        push!(nucleotides, "*" * uppercase(chain02))
        push!(sizes, length(chain02) + 1)
        push!(output_cards, html_div(
            children=[
                html_div("Second sequence type", className="title"),
                html_br(),
                html_div(getIdentity(uppercase(chain02)), className="content"),
            ], className="mini-card")
        )
    end

    kind_sequences = ""
    if getIdentity(uppercase(chain01)) == "PROTEIN" || getIdentity(uppercase(chain02)) == "PROTEIN" 
        kind_sequences = "PROTEIN"
    elseif getIdentity(uppercase(chain01)) == "RNA" || getIdentity(uppercase(chain02)) == "RNA"
        kind_sequences = "RNA"
    else
        kind_sequences = "DNA"
    end

    push!(output_cards, html_div(
            children=[
                html_div("SEQUENCES ALIGNED AS", className="title"),
                html_br(),
                html_div(kind_sequences, className="content"),
            ], className="mini-card")
        )

    global_matrix_scores, global_matrix_alignments = getMatricesScoreAlignments(nucleotides, sizes, kind_sequences)

    mm = findmax(sum(global_matrix_scores, dims=2))
    alignment_choose = string(mm[2][1]) * "th sequence choosen => " * nucleotides[mm[2][1]][begin+1:end] * " => with score " * string(mm[1])

    max_length_final_aligment, final_aligments, pairs_alignments_tags = getDifferentPairsAlignments(global_matrix_alignments, nucleotides, mm[2][1])

    # creating first and base alignment
    first_align = getFirstAlign(final_aligments, max_length_final_aligment)

    # creating msa
    msa = getMSA(final_aligments, first_align)

    # show(stdout, "text/plain", msa)
    return (output_cards, alignment_choose, pairs_alignments_tags, html_table(children=msa))
end # callback!

# Callback for secondary structure
callback!(
    app,
    # Output("out-plot", "children"),
    Output("max_value", "children"),
    Output("cytoscape-elements-callbacks", "children"),
    Input("button_secondary_structure", "n_clicks"),
    State("input_06", "value"),
    prevent_initial_call=true,) do n_clicks, chain01

    sizes = [length(chain01) + 1]
    nucleotides = ["*" * uppercase(chain01)]

    n = copy(nucleotides)
    n[1] = n[1][begin+1:end]
    matrix_scores_sa, matrix_directions_sa = create_matrix_sa(sizes[1] - 1)
    find_min = score_and_directions_sa(matrix_scores_sa, matrix_directions_sa, n, sizes[1] - 1)
    # println("-")
    # show(stdout, "text/plain", matrix_scores_sa)
    # println("-")
    # show(stdout, "text/plain", matrix_directions_sa)
    dynamic_plot = generate_secundary_structure(matrix_directions_sa, find_min, n)

    return (
        "Minimun value is $(find_min[1])",
        dynamic_plot
    )
end # callback!

# WEBAPP running
run_server(app, "0.0.0.0", debug=true)
