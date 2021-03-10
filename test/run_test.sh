if [ ! -d data ]; then
  echo "retrieving data"
  curl -OL https://www.dropbox.com/s/3r6ru1ww6evkw4w/grasp_data.tar.bz2
  tar xvjf grasp_data.tar.bz2
  rm -f grasp_data.tar.bz2
fi

echo "run example-1 (grasp)"
if [ ! -d output ]; then
  mkdir output
fi

DATA=data/samples

../bin/grasp -m 200 -n 1000 -o output ${DATA}/test.sdf ${DATA}/test.fasta ${DATA}/test.interaction 

echo "run example-2 (finder)"
if [ ! -d finder_out ]; then
  mkdir finder_out
fi

../bin/finder -o ./finder_out output/pair.patterns output/atom.dict ${DATA}/test.sdf ${DATA}/test.fasta ${DATA}/test.interaction

echo "run example-3 (evaluator)"
../bin/evaluator -o evaluate.dat output/pair.patterns output/atom.dict ${DATA}/random.sdf ${DATA}/random.fasta ${DATA}/random.pairs finder_out/feat.table

echo "run example-4 (predictor)"
../bin/predictor -n 100 -o predict.dat output/pair.patterns output/atom.dict ${DATA}/random.sdf ${DATA}/random.fasta finder_out/feat.table


