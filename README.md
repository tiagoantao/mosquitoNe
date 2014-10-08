mosquitoNe
==========

mosquitoNe

Requires Python 2 (not 3), matplotlib, NeEstimator2, SciPy

Very basic instructions
-----------------------


1. Run src/sim.py with a configuration file from conf, e.g.
python src/sim.py simple500

2. Concatenate genepop files
python src/concatenateGP.py simple500

3. Do stats:
python src/stats.py simple500

4. Figure plotting:
decline
python scripts/figDecline.py 60 20
python scripts/figDecline.py 60 50
wave plot
python scripts/plotWave.py season-500-700 60 20
python scripts/plotWave.py season-500-700 60 50


MLNE
----

1. Run
python src/mlne_stats.py simple500
(and all the other models: season-500-700 and decline-1000-100)
