# Bioinformatics Application

Bioinformatic application is a simple solution for: (1) global alignment, (2) local alignment, (3) star alignment, and (4) graphic generator of secondary structure.

It was created using a really nice programming language [julia](https://julialang.org/), [Dash](https://dash.plotly.com/julia/introduction) as framework, and docker for deployment.
Things that you will can see:

- [x] Quantity of characters of input strings
- [x] Identify is input string is DNA, RNA, or protein

## Main functions

### Alignments

- [x] Global Alignment
- [x] Local Alignment
- [x] Star Alignment

### Graphics

- [x] Secondary structure

### Some characteristics

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

- [ ] Add BLOSUM 50 80 90
- [ ] Scalable MSA in render frontend
- [ ] Download graphic of secondary structure
