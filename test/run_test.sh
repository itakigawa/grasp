echo "run example-1 (grasp)"
mkdir output
../bin/grasp -m 200 -n 1000 -o output ../data/test.sdf ../data/test.fasta ../data/test.interaction 

echo "run example-2 (finder)"
mkdir finder_out
../bin/finder -o ./finder_out output/pair.patterns output/atom.dict ../data/test.sdf ../data/test.fasta ../data/test.interaction

echo "run example-3 (evaluator)"
../bin/evaluator -o evaluate.dat output/pair.patterns output/atom.dict ../data/random.sdf ../data/random.fasta ../data/random.pairs finder_out/feat.table

echo "run example-4 (predictor)"
../bin/predictor -n 100 -o predict.dat output/pair.patterns output/atom.dict ../data/random.sdf ../data/random.fasta finder_out/feat.table


