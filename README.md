# Listing all frequent subgraph-substring pairs in graph-string pairs.

by Ichigaku Takigawa, Koji Tsuda, and Hiroshi Mamitsuka

## Paper

Takigawa I, Tsuda K, Mamitsuka H, 
Mining significant substructure pairs for interpreting polypharmacology in drug-target network. *PLoS One*, 2011 Feb 23;6(2):e16999.  doi: [10.1371/journal.pone.0016999](https://doi.org/10.1371/journal.pone.0016999).

[https://www.bic.kyoto-u.ac.jp/pathway/grasp/](https://www.bic.kyoto-u.ac.jp/pathway/grasp/)

## Compile and Install

Linux/Mac

```
$ git clone https://github.com/itakigawa/grasp.git
$ cd grasp/src
$ make; make install
```

Docker

```
$ git clone https://github.com/itakigawa/grasp.git
$ cd grasp
$ docker build -t grasp_dev .
$ docker run -it --rm  -v "$PWD":/workspace grasp_dev bash
bash# cd /workspace/src/
bash# make
bash# make install
```

## Input Format

Molecular graphs need to be encoded in [DFS code](HowToReadDFSCode.txt). See the gSpan paper for the detail. 

- Yan X, Han J, gSpan: Graph-Based Substructure Pattern Mining. Proc. 2002 of Int. Conf. on Data Mining (ICDM'02).

## Example

Data Retrieval

```
$ curl -OL https://www.dropbox.com/s/3r6ru1ww6evkw4w/grasp_data.tar.bz2
$ tar xvjf grasp_data.tar.bz2
```


### Example 1

Listing 1000 most significant substructure pairs shared by at least 200 interactions.

**Command:**

```
cd grasp/test
mkdir output
../bin/grasp -m 200 -n 1000 -o output test.sdf test.fasta test.interaction 
```

**Input:**   
- compound-protein pair: test.interaction
- compounds: test.sdf (the first line must include the ID)
- proteins: test.fasta (the first line must include the ID)

**Output:**
- outout/pair.patterns   ... 10000 most significant substructure pairs
- outout/feat.instances  ... i-th line shows compound-protein pairs containing i-th substructure pairs
- outout/atom.dict ... dictionary for DFS code

### Example 2

Profiling arbitrary compound-ligand pairs using obtained substructure pairs.

```
$ cd grasp/test
$ mkdir finder_out
$ ../bin/finder -o ./finder_out output/pair.patterns output/atom.dict test.sdf test.fasta test.interaction
```

*Output:*
- ./finder_out/feat.table, the i-th line of which represents the substructure pairs contained in the i-th compound-protein pair.

### Example 3

Nearest neighbor search for pairs explicitly defined.

```
$ ../bin/evaluator -o evaluate.dat output/pair.patterns output/atom.dict random.sdf random.fasta random.pairs finder_out/feat.table
```

100 highest scoring nearest neighbor search for all possible combinations of pairs.

```
$ ../bin/predictor -n 100 -o predict.dat output/pair.patterns output/atom.dict random.sdf random.fasta finder_out/feat.table
```


## Paper's Case

To obtain fragment-pair patterns in the paper

```
$ cd grasp
$ ./bin/grasp -m 561 -n 20000 -o grasp_out data/drugbank_sm/drug.sdf data/drugbank_sm/target.fasta data/drugbank_sm/interaction.txt > grasp_out/log | cat
$ ruby bin/feat_reductor.rb grasp_out/feat.instances grasp_out/pair.patterns 10000 > nonredundant_10000.patterns
```

To generate feat.table

```
$ time ./bin/finder -o ./finder_out nonredundant_10000.patterns grasp_out/atom.dict data/drugbank_sm/drug.sdf data/drugbank_sm/target.fasta data/drugbank_sm/interaction.txt > finder_out/log
```
