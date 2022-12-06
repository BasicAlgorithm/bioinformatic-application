module UTILITIES

using Dash
using DashHtmlComponents
using DashCoreComponents

export dropdown_initial_options
export dropdown_q_msa
export functions_tags
export input_msa

dropdown_initial_options = [
  Dict("label" => "Global and local alignment", "value" => "NYC"),
  Dict("label" => "Star alignment - MSA", "value" => "MTL"),
  Dict("label" => "Secondary graphic structure", "value" => "SF"),
]

dropdown_q_msa = [
  Dict("label" => "2", "value" => 2),
  Dict("label" => "3", "value" => 3),
  Dict("label" => "4", "value" => 4),
  Dict("label" => "5", "value" => 5),
]

functions_tags = Dict([
  ("NYC", [
    html_div(
      children=[
        dcc_input(id="input_01", type="text", value="AGC", className="form__field"),
        html_label("Write your first sequence", className="form__label"),
      ], className="form__group field"),
    html_div(
      children=[
        dcc_input(id="input_02", type="text", value="AAAC", className="form__field"),
        html_label("Write your second sequence", className="form__label"),
      ], className="form__group field"),
    html_div(children=[
        html_p("Choose maximun alignments"),
        html_br(),
        html_br(),
        dcc_slider(
          id="slider-value",
          min=1,
          max=8,
          step=nothing,
          value=1,
          marks=Dict(
            1 => Dict("label" => "1"),
            2 => Dict("label" => "2"),
            3 => Dict("label" => "3"),
            4 => Dict("label" => "4"),
            5 => Dict("label" => "5"),
            6 => Dict("label" => "6"),
            7 => Dict("label" => "7"),
            8 => Dict("label" => "8"),
          )
        )
      ], className="slider"),
    html_div(
      children=[
        html_button(id="button_global_local_alig", children="Analize", className="button-89", n_clicks=0),
      ], className="div_button"),
    html_br(),
    html_br(),
    html_br(),
    html_div(
      html_h1("GLOBAL ALIGNMENT"), className="card-title"),
    html_div(
      children=[
        html_div(
          children=[
            html_div("First sequence lenght", className="title"),
            html_br(),
            html_div(id="lenght_01", className="content"),
          ], className="mini-card"),
        html_div(
          children=[
            html_div("Second sequence lenght", className="title"),
            html_br(),
            html_div(id="lenght_02", className="content"),
          ], className="mini-card"),
        html_div(
          children=[
            html_div("First sequence type", className="title"),
            html_br(),
            html_div(id="identified_01", className="content"),
          ], className="mini-card"),
        html_div(
          children=[
            html_div("Second sequence type", className="title"),
            html_br(),
            html_div(id="identified_02", className="content"),
          ], className="mini-card"),
      ], className="mini-cards"),
    html_div(
      children=[
        html_div("Alignments", className="title"),
        html_br(),
        html_div(id="alignments", children=[]),
      ], className="card"),
    html_div(
      children=[
        html_div("Final score", className="title"),
        html_br(),
        html_div(id="score", className="content"),
      ], className="card"),
    html_div(
      html_h1("LOCAL ALIGNMENT"), className="card-title"),
    html_div(
      children=[
        html_div(id="la_aligments"),
      ], className="card"),
  ]),
  ("MTL", [
    html_div([
        dcc_dropdown(id="q_sequence_msa", options=dropdown_q_msa, placeholder="choose quantity of sequences", value="")
      ], className="dropdown"),
    html_div(id="msa_with_multi_sequence"),]),
  ("SF", [
    html_div(
      children=[
        dcc_input(id="input_06", type="text", value="GGAACUAUC", className="form__field"),
        html_label("Write your sequence", className="form__label"),
      ], className="form__group field"),
    html_div(
      children=[
        html_button(id="button_secondary_structure", children="create graphic structure", className="button-89", n_clicks=0),
      ], className="div_button"),
    html_div(
      html_h1("SECONDARY STRUCTURE"), className="card-title"),
    html_div(id="max_value", className="card"),
    html_div(id="cytoscape-elements-callbacks", className="card_cytoscape"),
  ]),
])

input_msa = [
  html_div(
    children=[
      dcc_input(id="input_03", type="text", value="ACCCGTU", className="form__field"),
      html_label("Write your first sequence", className="form__label"),
    ], className="form__group field"),
  html_div(
    children=[
      dcc_input(id="input_04", type="text", value="AAAC", className="form__field"),
      html_label("Write your second sequence", className="form__label"),
    ], className="form__group field"),
  html_div(
    children=[
      dcc_input(id="input_05", type="text", value="GCAA", className="form__field"),
      html_label("Write your third sequence", className="form__label"),
    ], className="form__group field"),
  html_div(
    children=[
      html_button(id="button_star_alignment", children="Align sequences", className="button-89", n_clicks=0),
    ], className="div_button"),
  html_br(),
  html_br(),
  html_br(),
  html_div(
    children=[
      html_div(
        children=[
          html_div("First sequence type", className="title"),
          html_br(),
          html_div(id="identified_03", className="content"),
        ], className="mini-card"),
      html_div(
        children=[
          html_div("Second sequence type", className="title"),
          html_br(),
          html_div(id="identified_04", className="content"),
        ], className="mini-card"),
      html_div(
        children=[
          html_div("Third sequence type", className="title"),
          html_br(),
          html_div(id="identified_05", className="content"),
        ], className="mini-card"),
      html_div(
        children=[
          html_div("Aligned as", className="title"),
          html_br(),
          html_div(id="identified_06", className="content"),
        ], className="mini-card"),
    ], className="mini-cards"),
  html_div(
    html_h1("Alignment choosen"), className="card-title"),
  html_div(
    children=[
      html_div(id="sa_align_choosen"),
    ], className="card"),
  html_div(
    html_h1("Pair alignments individually"), className="card-title"),
  html_div(
    children=[
      html_div(id="sa_pairs_aligments"),
    ], className="card"),
  html_div(
    html_h1("Multiple alignment MSA"), className="card-title"),
  html_div(
    children=[
      html_div(id="sa_multiple_aligments"),
    ], className="card"),
]

end