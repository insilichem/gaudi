exec("""
import sys
sys.path[:0] = ['/home/jr/dev/gaudi']
import gaudi
import deap
from deap import creator, tools
cfg = gaudi.parse.Settings('/home/jr/dev/modulartest.yaml')
gaudi.plugin.import_plugins(*cfg.genes)
gaudi.plugin.import_plugins(*cfg.objectives)
toolbox = deap.base.Toolbox()
deap.creator.create("Fitness", gaudi.base.Fitness,objectivelist=cfg.objectives,weights=cfg.weights)
deap.creator.create("Individual", gaudi.base.Individual,fitness=deap.creator.Fitness)
toolbox.register("call", (lambda fn, *args, **kwargs: fn(*args, **kwargs)))
toolbox.register("individual", toolbox.call, deap.creator.Individual, cfg.genes, cfg=cfg)
ind = toolbox.individual()
ind.express()
""")