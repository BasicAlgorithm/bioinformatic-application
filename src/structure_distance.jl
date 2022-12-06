# method of secondary Structure
# https://www.frontiersin.org/articles/10.3389/fgene.2019.00467/full

using SparseArrays: SparseMatrixCSC, sparse
using ArnoldiMethod: SR
using Base: OneTo
using LinearAlgebra: eigen

function random_layout(g)
  rand(nv(g)), rand(nv(g))
end

function circular_layout(g)
  if nv(g) == 1
    return [0.0], [0.0]
  else
    # Discard the extra angle since it matches 0 radians.
    θ = range(0, stop=2pi, length=nv(g) + 1)[1:end-1]
    return cos.(θ), sin.(θ)
  end
end

function spring_layout(g::AbstractGraph,
  locs_x=2 * rand(nv(g)) .- 1.0,
  locs_y=2 * rand(nv(g)) .- 1.0;
  C=2.0,
  MAXITER=100,
  INITTEMP=2.0)

  nvg = nv(g)
  adj_matrix = adjacency_matrix(g)

  # The optimal distance bewteen vertices
  k = C * sqrt(4.0 / nvg)
  k² = k * k

  # Store forces and apply at end of iteration all at once
  force_x = zeros(nvg)
  force_y = zeros(nvg)

  # Iterate MAXITER times
  @inbounds for iter = 1:MAXITER
    # Calculate forces
    for i = 1:nvg
      force_vec_x = 0.0
      force_vec_y = 0.0
      for j = 1:nvg
        i == j && continue
        d_x = locs_x[j] - locs_x[i]
        d_y = locs_y[j] - locs_y[i]
        dist² = (d_x * d_x) + (d_y * d_y)
        dist = sqrt(dist²)

        if !(iszero(adj_matrix[i, j]) && iszero(adj_matrix[j, i]))
          # Attractive + repulsive force
          # F_d = dist² / k - k² / dist # original FR algorithm
          F_d = dist / k - k² / dist²
        else
          # Just repulsive
          # F_d = -k² / dist  # original FR algorithm
          F_d = -k² / dist²
        end
        force_vec_x += F_d * d_x
        force_vec_y += F_d * d_y
      end
      force_x[i] = force_vec_x
      force_y[i] = force_vec_y
    end
    # Cool down
    temp = INITTEMP / iter
    # Now apply them, but limit to temperature
    for i = 1:nvg
      fx = force_x[i]
      fy = force_y[i]
      force_mag = sqrt((fx * fx) + (fy * fy))
      scale = min(force_mag, temp) / force_mag
      locs_x[i] += force_x[i] * scale
      locs_y[i] += force_y[i] * scale
    end
  end

  # Scale to unit square
  min_x, max_x = minimum(locs_x), maximum(locs_x)
  min_y, max_y = minimum(locs_y), maximum(locs_y)
  function scaler(z, a, b)
    2.0 * ((z - a) / (b - a)) - 1.0
  end
  map!(z -> scaler(z, min_x, max_x), locs_x, locs_x)
  map!(z -> scaler(z, min_y, max_y), locs_y, locs_y)

  return locs_x, locs_y
end

using Random: MersenneTwister

function spring_layout(g::AbstractGraph, seed::Integer, kws...)
  rng = MersenneTwister(seed)
  spring_layout(g, 2 .* rand(rng, nv(g)) .- 1.0, 2 .* rand(rng, nv(g)) .- 1.0; kws...)
end

"""
This function is copy from [IainNZ](https://github.com/IainNZ)'s [GraphLayout.jl](https://github.com/IainNZ/GraphLayout.jl)
"""

function shell_layout(g, nlist::Union{Nothing,Vector{Vector{Int}}}=nothing)
  if nv(g) == 1
    return [0.0], [0.0]
  end
  if nlist === nothing
    nlist = [collect(1:nv(g))]
  end
  radius = 0.0
  if length(nlist[1]) > 1
    radius = 1.0
  end
  locs_x = Float64[]
  locs_y = Float64[]
  for nodes in nlist
    # Discard the extra angle since it matches 0 radians.
    θ = range(0, stop=2pi, length=length(nodes) + 1)[1:end-1]
    append!(locs_x, radius * cos.(θ))
    append!(locs_y, radius * sin.(θ))
    radius += 1.0
  end
  return locs_x, locs_y
end

function spectral_layout(g::AbstractGraph, weight=nothing)
  if nv(g) == 1
    return [0.0], [0.0]
  elseif nv(g) == 2
    return [0.0, 1.0], [0.0, 0.0]
  end

  if weight === nothing
    weight = ones(ne(g))
  end
  if nv(g) > 500
    A = sparse(Int[src(e) for e in edges(g)],
      Int[dst(e) for e in edges(g)],
      weight, nv(g), nv(g))
    if is_directed(g)
      A = A + transpose(A)
    end
    return _spectral(A)
  else
    L = laplacian_matrix(g)
    return _spectral(Matrix(L))
  end
end

function _spectral(L::Matrix)
  eigenvalues, eigenvectors = eigen(L)
  index = sortperm(eigenvalues)[2:3]
  return eigenvectors[:, index[1]], eigenvectors[:, index[2]]
end

function _spectral(A::SparseMatrixCSC)
  data = vec(sum(A, dims=1))
  D = sparse(Base.OneTo(length(data)), Base.OneTo(length(data)), data)
  L = D - A
  eigenvalues, eigenvectors = Graphs.LinAlg.eigs(L, nev=3, which=SR())
  index = sortperm(real(eigenvalues))[2:3]
  return real(eigenvectors[:, index[1]]), real(eigenvectors[:, index[2]])
end

function gplot(g::AbstractGraph{T},
  locs_x_in::Vector{R1}, locs_y_in::Vector{R2};
  nodelabel = nothing,
  nodelabelc = colorant"black",
  nodelabelsize = 1.0,
  NODELABELSIZE = 4.0,
  nodelabeldist = 0.0,
  nodelabelangleoffset = π / 4.0,
  edgelabel = [],
  edgelabelc = colorant"black",
  edgelabelsize = 1.0,
  EDGELABELSIZE = 4.0,
  edgestrokec = colorant"lightgray",
  edgelinewidth = 1.0,
  EDGELINEWIDTH = 3.0 / sqrt(nv(g)),
  edgelabeldistx = 0.0,
  edgelabeldisty = 0.0,
  nodesize = 1.0,
  NODESIZE = 0.25 / sqrt(nv(g)),
  nodefillc = colorant"turquoise",
  nodestrokec = nothing,
  nodestrokelw = 0.0,
  arrowlengthfrac = is_directed(g) ? 0.1 : 0.0,
  arrowangleoffset = π / 9,
  linetype = "straight",
  outangle = π / 5) where {T <:Integer, R1 <: Real, R2 <: Real}

  length(locs_x_in) != length(locs_y_in) && error("Vectors must be same length")
  N = nv(g)
  NE = ne(g)
  if nodelabel != nothing && length(nodelabel) != N
      error("Must have one label per node (or none)")
  end
  if !isempty(edgelabel) && length(edgelabel) != NE
      error("Must have one label per edge (or none)")
  end

  locs_x = Float64.(locs_x_in)
  locs_y = Float64.(locs_y_in)

  # Scale to unit square
  min_x, max_x = extrema(locs_x)
  min_y, max_y = extrema(locs_y)
  function scaler(z, a, b)
      if (a - b) == 0.0
          return 0.5
      else
          return 2.0 * ((z - a) / (b - a)) - 1.0
      end
  end
  map!(z -> scaler(z, min_x, max_x), locs_x, locs_x)
  map!(z -> scaler(z, min_y, max_y), locs_y, locs_y)

  # Determine sizes
  #NODESIZE    = 0.25/sqrt(N)
  #LINEWIDTH   = 3.0/sqrt(N)

  max_nodesize = NODESIZE / maximum(nodesize)
  nodesize *= max_nodesize
  max_edgelinewidth = EDGELINEWIDTH / maximum(edgelinewidth)
  edgelinewidth *= max_edgelinewidth
  max_edgelabelsize = EDGELABELSIZE / maximum(edgelabelsize)
  edgelabelsize *= max_edgelabelsize
  max_nodelabelsize = NODELABELSIZE / maximum(nodelabelsize)
  nodelabelsize *= max_nodelabelsize
  max_nodestrokelw = maximum(nodestrokelw)
  if max_nodestrokelw > 0.0
      max_nodestrokelw = EDGELINEWIDTH / max_nodestrokelw
      nodestrokelw *= max_nodestrokelw
  end

  # Create nodes
  nodecircle = fill(0.4Compose.w, length(locs_x))
  if isa(nodesize, Real)
            for i = 1:length(locs_x)
                  nodecircle[i] *= nodesize
            end
    else
          for i = 1:length(locs_x)
                  nodecircle[i] *= nodesize[i]
            end
    end
  nodes = circle(locs_x, locs_y, nodecircle)

  # Create node labels if provided
  texts = nothing
  if nodelabel != nothing
      text_locs_x = deepcopy(locs_x)
      text_locs_y = deepcopy(locs_y)
      texts = text(text_locs_x .+ nodesize .* (nodelabeldist * cos(nodelabelangleoffset)),
                   text_locs_y .- nodesize .* (nodelabeldist * sin(nodelabelangleoffset)),
                   map(string, nodelabel), [hcenter], [vcenter])
  end
  # Create edge labels if provided
  edgetexts = nothing
  if !isempty(edgelabel)
      edge_locs_x = zeros(R, NE)
      edge_locs_y = zeros(R, NE)
      for (e_idx, e) in enumerate(edges(g))
          i = src(e)
          j = dst(e)
          mid_x = (locs_x[i]+locs_x[j]) / 2.0
          mid_y = (locs_y[i]+locs_y[j]) / 2.0
          edge_locs_x[e_idx] = (is_directed(g) ? (mid_x+locs_x[j]) / 2.0 : mid_x) + edgelabeldistx * NODESIZE
          edge_locs_y[e_idx] = (is_directed(g) ? (mid_y+locs_y[j]) / 2.0 : mid_y) + edgelabeldisty * NODESIZE

      end
      edgetexts = text(edge_locs_x, edge_locs_y, map(string, edgelabel), [hcenter], [vcenter])
  end

  # Create lines and arrow heads
  lines, arrows = nothing, nothing
  if linetype == "curve"
      if arrowlengthfrac > 0.0
          curves_cord, arrows_cord = graphcurve(g, locs_x, locs_y, nodesize, arrowlengthfrac, arrowangleoffset, outangle)
          lines = curve(curves_cord[:,1], curves_cord[:,2], curves_cord[:,3], curves_cord[:,4])
          arrows = line(arrows_cord)
      else
          curves_cord = graphcurve(g, locs_x, locs_y, nodesize, outangle)
          lines = curve(curves_cord[:,1], curves_cord[:,2], curves_cord[:,3], curves_cord[:,4])
      end
  else
      if arrowlengthfrac > 0.0
          lines_cord, arrows_cord = graphline(g, locs_x, locs_y, nodesize, arrowlengthfrac, arrowangleoffset)
          lines = line(lines_cord)
          arrows = line(arrows_cord)
      else
          lines_cord = graphline(g, locs_x, locs_y, nodesize)
          lines = line(lines_cord)
      end
  end

  compose(context(units=UnitBox(-1.2, -1.2, +2.4, +2.4)),
          compose(context(), texts, fill(nodelabelc), stroke(nothing), fontsize(nodelabelsize)),
          compose(context(), nodes, fill(nodefillc), stroke(nodestrokec), linewidth(nodestrokelw)),
          compose(context(), edgetexts, fill(edgelabelc), stroke(nothing), fontsize(edgelabelsize)),
          compose(context(), arrows, stroke(edgestrokec), linewidth(edgelinewidth)),
          compose(context(), lines, stroke(edgestrokec), fill(nothing), linewidth(edgelinewidth)))
end

function gplot(g; layout::Function=spring_layout, keyargs...)
  gplot(g, layout(g)...; keyargs...)
end

# take from [Gadfly.jl](https://github.com/dcjones/Gadfly.jl)
function open_file(filename)
  if Sys.isapple() #apple
      run(`open $(filename)`)
  elseif Sys.islinux() || Sys.isbsd() #linux
      run(`xdg-open $(filename)`)
  elseif Sys.iswindows() #windows
      run(`$(ENV["COMSPEC"]) /c start $(filename)`)
  else
      @warn("Showing plots is not supported on OS $(string(Sys.KERNEL))")
  end
end

# taken from [Gadfly.jl](https://github.com/dcjones/Gadfly.jl)
function gplothtml(args...; keyargs...)
  filename = string(tempname(), ".html")
  output = open(filename, "w")

  plot_output = IOBuffer()
  draw(SVGJS(plot_output, Compose.default_graphic_width,
             Compose.default_graphic_width, false), gplot(args...; keyargs...))
  plotsvg = String(take!(plot_output))

  write(output,
      """
      <!DOCTYPE html>
      <html>
        <head>
          <title>GraphPlot Plot</title>
          <meta charset="utf-8">
        </head>
          <body>
          <script charset="utf-8">
              $(read(Compose.snapsvgjs, String))
          </script>
          <script charset="utf-8">
              $(read(gadflyjs, String))
          </script>
          $(plotsvg)
        </body>
      </html>
      """)
  close(output)
  open_file(filename)
end