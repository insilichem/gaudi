#!/usr/bin/python

# MOF3D
# Multi Objective Force Field Free Docking
# A docking module for UCSF Chimera
# Jaime RGP <https://bitbucket.org/jrgp> @ UAB, 2014

class Settings(object):
	''' Simple parser for INI settings file. Allows both dict and object syntaxes '''
	def __init__(self, path, asDict=False):
		self._path = path
		self._parse(bool(asDict))
	
	class Param(object):
		def __init__(self, *d):
			for d_ in d:
				self.__dict__.update(d_)

	def _parse(self, asDict):
		with open(self._path, 'r') as config:
			s = None
			parsed = {}
			for line in config.readlines():
				line = line.strip()
				if line == '' or line.startswith('#'):
					continue
				elif line.startswith('[') and line.endswith(']'):
					s = line.strip('[]')
					if s in parsed and isinstance(parsed[s], dict):
						parsed[s] = [parsed[s]]
					elif s in parsed and isinstance(parsed[s], list):
						parsed[s].append({})
					else:
						parsed.setdefault(s, dict())
				elif '=' in line:
					line = line.split('#',1)[0]
					k, v = [ _.strip() for _ in line.split('=',1) ]
					if ', ' in v:
						v = [ self._num(x) for x in v.split(', ')]
					else:
						v = self._num(v)

					if isinstance(parsed[s], dict):
						parsed[s][k] = v
					elif isinstance(parsed[s], list):
						parsed[s][-1][k] = v	
		if asDict:
			self.parsed = parsed
		self._to_object(parsed)
		del parsed

	def _to_object(self, d):
		for s, c in d.items():
			if isinstance(c, dict):
				setattr(self, s, self.Param(c))
			if isinstance(c, list):
				setattr(self, s, self.Param({'number': len(c)}))
				for i, l in enumerate(c):
					l['priority'] = i
					name = l['feature']
					setattr(getattr(self, s), name, self.Param(l))
	def _num(self, n):
		try:
			if '.' in n or 'e' in n or 'E' in n:
				return float(n)
			else:
				return int(n)
		except ValueError:
			return n

#############################
# def settings(path):
# 	with open(path, 'r') as config:
# 		s = None
# 		parsed = {}
# 		for line in config.readlines():
# 			line = line.strip()
# 			if line == '' or line.startswith('#'):
# 				continue
# 			elif line.startswith('[') and line.endswith(']'):
# 				s = line.strip('[]')
# 				if s in parsed and isinstance(parsed[s], dict):
# 					parsed[s] = [parsed[s]]
# 				elif s in parsed and isinstance(parsed[s], list):
# 					parsed[s].append({})
# 				else:
# 					parsed.setdefault(s, dict())
# 			elif '=' in line:
# 				k, v = [ _.strip() for _ in line.split('=',1) ]
# 				if ', ' in v:
# 					v = [ _num(x) for x in v.split(', ')]
# 				else:
# 					v = _num(v)

# 				if isinstance(parsed[s], dict):
# 					parsed[s][k] = v
# 				elif isinstance(parsed[s], list):
# 					parsed[s][-1][k] = v
# 	config.close()
# 	return parsed				

def _test_rebuild(cfg):
	for s, c in cfg.items():
		if isinstance(c, dict):
			print '\n['+s+']'
			for k, v in c.items():
				print k, "=", v
		elif isinstance(c, list):
			for i, l in enumerate(c):
				print '\n['+str(s)+' '+str(i)+']'
				for k, v in l.items():
					print k, "=", v

if __name__ == '__main__':
	#cfg = settings('/home/jr/x/hyde/mof3d.ini')
	settings = Settings('/home/jr/x/hyde/mof3d.ini')
	print (settings.rotamers.residues)

