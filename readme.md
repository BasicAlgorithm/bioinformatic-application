# Bioinformatics Application

Bioinformatic application is a simple solution for: (1) global alignment, (2) local alignment, (3) Multiple Sequence Alignment MSA using star alignment, and (4) graphic generator of secondary structure.

It was created using a really nice programming language [julia](https://julialang.org/), [Dash](https://dash.plotly.com/julia/introduction) as framework, and docker for deployment.

## Main functions

### Alignments

- [x] Global Alignment | Neddleman-Wunch algorithm
- [x] Local Alignment | Smith-Waterman algorithm
- [x] Multiple Sequence Alignment MSA | Gusfield Star alignment

### Graphics

- [x] Secondary structure prediction of RNA | Watson-Crick algorithm

### Some characteristics

- [x] Length of input strings.
- [x] Identify if input string is DNA, RNA, or protein.
- [x] BLOSUM62 and Standard ADN ARN table score
- [x] Alignments have colors for better differentiation
- [x] Length and type of alignments
- [x] Score of alignments
- [x] Dynamic secondary structure

## Run the project

You have to have installed docker and docker-compose.

```bash
docker-compose build
docker-compose up
```

Now, open localhost:3000 or [click here](http://localhost:3000/). If this port is busy, you can change the port on docker-compose.yaml file.

## TODO

- [ ] Add BLOSUM 50 80 90.
- [ ] Scalable MSA in render frontend.
- [ ] Download graphic of secondary structure.
- [ ] Restrict secondary structure for only ARN, otherwise show error.
- [ ] secondary structure: new unions with different color.
