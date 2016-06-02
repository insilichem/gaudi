from gaudi.base import Individual, Environment
from gaudi.parse import Settings
import logging
l = logging.getLogger()
l.addHandler(logging.StreamHandler())

cfg = Settings(path="/path/to/some/test.gaudi-input")

ind = Individual(cfg=cfg)
ind.express()
ind.unexpress()

gene = ind.genes['SomeGeneName']
gene.express()
gene.unexpress()

env = Environment(cfg)
env.evaluate(ind)

obj = env.objectives['SomeObjectiveName']
obj.evaluate(ind)