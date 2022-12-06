FROM julia:1.8.3-buster

ARG APP_HOME=/dash

WORKDIR ${APP_HOME}

RUN julia -e "using Pkg; Pkg.add.(\"Dash\");"
RUN julia -e "using Pkg; Pkg.add.(\"DashHtmlComponents\");"
RUN julia -e "using Pkg; Pkg.add.(\"DashCoreComponents\");"
RUN julia -e "using Pkg; Pkg.add.(\"Compose\");"
RUN julia -e "using Pkg; Pkg.add.(\"Colors\");"
RUN julia -e "using Pkg; Pkg.add.(\"Cairo\");"
RUN julia -e "using Pkg; Pkg.add.(\"Fontconfig\");"
RUN julia -e "using Pkg; Pkg.add.(\"Graphs\");"
RUN julia -e "using Pkg; Pkg.add.(\"GraphPlot\");"
RUN julia -e "using Pkg; Pkg.add.(\"DashCytoscape\");"

RUN julia -e "using Pkg; Pkg.update;"

COPY . .

EXPOSE 8050

CMD ["julia", "src/app.jl"]