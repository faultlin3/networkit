#
# file: profiling.py
# author: Mark Erb
#

from networkit import *
import networkit as kit

import os as os
import sys, traceback

from . import multiprocessing
from . import stat
from . import plot

from IPython.core.display import *
import collections
import math
import fnmatch


def readfile(filename, removeWS=False):
	""" private helper function for file-loading """
	with open(__file__[:__file__.rfind("/")+1] + filename, "r") as file:
		result = file.read()
		if removeWS:
			result = " ".join(result.split())
		return result


try:
	__IPYTHON__

	def _initHeader(tag, type, data):
		""" private helper function for ipython/jupyther-notebook hack: create content of extended header """
		result = """
			{
				var element = document.getElementById('NetworKit_""" + tag + """');
				if (element) {
					element.parentNode.removeChild(element);
				}
				element = document.createElement('""" + tag + """');
				element.type = 'text/""" + type + """';
				element.innerHTML = '""" + data + """';
				element.setAttribute('id', 'NetworKit_""" + tag + """');
				document.head.appendChild(element);
			}
		"""
		return result


	def _initOverlay(name, data):
		""" private helper function for ipython/jupyther-notebook hack: create content of overlay """
		result = """
			{
				var element = document.getElementById('NetworKit_""" + name + """');
				if (element) {
					element.parentNode.removeChild(element);
				}
				element = document.createElement('div');
				element.innerHTML = '<div id="NetworKit_""" + name + """_Toolbar_Top"><div class="button icon-close" id="NetworKit_""" + name + """_Close" /></div>""" + data + """';
				element.setAttribute('id', 'NetworKit_""" + name + """');
				document.body.appendChild(element);
				document.getElementById('NetworKit_""" + name + """_Close').onclick = function (e) {
					document.getElementById('NetworKit_""" + name + """').style.display = 'none';
				}
			}
		"""
		return result

	# load headers/scripts
	display_html(
		HTML("""
			<script type="text/javascript">
			<!--
				""" + _initHeader("script", "javascript", readfile("html/profiling.js", True))  + """
				""" + _initHeader("style",  "css",        readfile("html/profiling.css", True)) + """
				""" + _initOverlay("Overlay", readfile("html/overlay.html", True)) + """
			-->
			</script>
		""")
	)
except NameError as e:
	pass


class Config:
	""" object to control profiling.Profile behaviour """
	
	def __init__(self):
		""" constructor: all options off """
		self.__options_Properties = {
			"Diameter": False
		}
		self.__options_Measures = {
			"Centrality.Degree": False,
			"Centrality.CoreDecomposition": False,
			"Centrality.ClusteringCoefficient": False,
			"Centrality.PageRank": False,
			"Centrality.KPath": False,
			"Centrality.Katz": False,
			"Centrality.Betweenness": False,
			"Centrality.Closeness": False,
			"Partition.Communities": False,
			"Partition.ConnectedComponents": False,
			"Partition.CoreDecomposition": False
		}
		self.__options_MeasureCorrelations = {
			"Pearson": False,
			"Spearman": False,
			"Fechner": False
		}

	@classmethod
	def createConfig(cls, preset="default"):
		""" create a config object by selection a named preset of predefined options (complete/minimal/default) """
		result = Config()

		if preset == "complete":
			result.setProperty("Diameter")
			result.setMeasure("Centrality.Degree"),
			result.setMeasure("Centrality.CoreDecomposition")
			result.setMeasure("Centrality.ClusteringCoefficient")
			result.setMeasure("Centrality.PageRank")
			result.setMeasure("Centrality.Katz")
			result.setMeasure("Centrality.Betweenness")
			result.setMeasure("Centrality.Closeness")
			result.setMeasure("Partition.Communities")
			result.setMeasure("Partition.ConnectedComponents")
			result.setMeasure("Partition.CoreDecomposition")
			result.setMeasureCorrelation("Spearman")
		
		elif preset == "minimal":
			result.setMeasure("Centrality.Degree")
			result.setMeasure("Partition.ConnectedComponents")
		
		elif preset == "default":
			result.setProperty("Diameter")
			result.setMeasure("Centrality.Degree")
			result.setMeasure("Centrality.ClusteringCoefficient")
			result.setMeasure("Centrality.PageRank")
			result.setMeasure("Centrality.Betweenness")
			result.setMeasure("Centrality.Katz")
			result.setMeasure("Centrality.CoreDecomposition")
			result.setMeasure("Partition.ConnectedComponents")
			result.setMeasure("Partition.Communities")
			result.setMeasure("Partition.CoreDecomposition")
			result.setMeasureCorrelation("Spearman")
		else:
			raise Error("no preset given")

		return result

		
	def setProperty(self, id, enabled=True):
		""" enable/disable given property option """
		if id in self.__options_Properties:
			self.__options_Properties[id] = enabled


	def getProperty(self, id):
		""" return state of given property option """
		return self.__options_Properties[id]


	def setMeasure(self, id, enabled=True):
		""" enable/disable given measure option """
		if id in self.__options_Measures:
			self.__options_Measures[id] = enabled


	def getMeasure(self, id):
		""" return state of given measure option """
		return self.__options_Measures[id]


	def setMeasureCorrelation(self, id, enabled=True):
		""" enable/disable given correlation (measure) option """
		if id in self.__options_MeasureCorrelations:
			self.__options_MeasureCorrelations[id] = enabled


	def getMeasureCorrelation(self, id):
		""" return state of given correlation (measure) option """
		return self.__options_MeasureCorrelations[id]



class Profile:
	""" automated network profiling class """
	
	__TOKEN = object()	# see __init__(): prevent this class from being instanced directly
	__pageCount = 0
	__verbose = False
	__verboseLevel = 0
	__verboseFilename = ""
	__parallel = multiprocessing.numberOfProcessors()

	
	def __init__(self, G, token=object()):
		""" constructor: don't call this method directly - use create instead """
		if token is not self.__TOKEN:
			raise ValueError("call create(G) to create an instance")
		self.__G = G
		self.__properties = {}
		self.__measures = collections.OrderedDict()
		self.__correlations = {}

		
	@classmethod
	def walk(cls, inputDir, filePattern="*", outputDir=None, config=Config(), outputType="HTML", style="light", color=(0, 0, 1), recursive=False, parallel=False, graphFormat=None):
		""" tests all files of a directory for the given conditions and generates a profile when matching
		
		Args:
			inputDir: the directory to search
			filePattern: specify accepted file names, e.g.: *.METIS.graph
			outputDir: directory to write the generated profiles
			config: object to control some aspects of the generation behaviour (Config)
			outputType: profile output format ("HTML", "LaTeX")
			style: style of generated output ("light")
			color: mainly used color of given style
			recursive: also search in subfolders for matching files
			parallel: run some additional parts of the generation in parallel (experimental)
			graphFormat: format of matching files (e.g.: Format.METIS)
		"""
		
		if outputDir is None:
			raise ValueError("output directory (parameter: outputDir) not specified")
		elif not os.path.isdir(outputDir):
			os.mkdir(outputDir)
		if graphFormat is None:
			raise ValueError("no graph format has been specified")
		
		for (dirpath, dirnames, filenames) in os.walk(inputDir):
			for filename in filenames:
				file = dirpath + "/" + filename
				if fnmatch.fnmatch(filename, filePattern):
					cls.__verbosePrint("\n[ " + file + " ]")
					try:
						G = kit.readGraph(file, graphFormat)
						try:
							pf = cls.create(
								G,
								config = config
							)
							cls.__verbosePrint("");
							pf.output(
								outputType = outputType,
								directory = outputDir,
								style = style,
								color = color,
								parallel = parallel
							)
						except Exception as e:
							cls.__verbosePrint("=> an error occured: {0} of type {1}".format(e, type(e)))
							cls.__verbosePrint(traceback.format_exc())
					except:
						cls.__verbosePrint("could not read {0}".format(file))
					cls.__verbosePrint("\n")
				else:
					cls.__verbosePrint("skipping {0} as it does not match filePattern".format(file))
			if not recursive:
				break
		print("Done")


	@classmethod
	def create(cls, G, config=Config()):
		""" creates a profile object
		
		Args:
			G: graph to profile
			config: object to control some aspects of the generation behaviour (Config)
		Returns:
			profile object
		"""

		result = cls(G, cls.__TOKEN)
		# TODO: use copy constructor instead
		result.__config = config

		kit.setNumberOfThreads(result.__parallel)

		def funcScores(instance):
			""" returns node scores """
			return instance.scores()

		def funcSizes(instance):
			""" returns partition subset sizes """
			return sorted(instance.getPartition().subsetSizes())

		if G.isDirected():
			classConnectedComponents = components.StronglyConnectedComponents
		else:
			classConnectedComponents = components.ConnectedComponents

		# internal unique name | category name | display name | 
		# compute correlation within same category | value function for measures | display name (axis) | class name of measure | parameter of constructor
		for parameter in [
			("Centrality.Degree",					"Node Centrality",	"Degree",
				True,	funcScores,	"Score",				centrality.DegreeCentrality, 			(G, )),
			("Centrality.CoreDecomposition",		"Node Centrality",	"k-Core Decomposition",
				True,	funcScores,	"Score",				centrality.CoreDecomposition, 			(G, )),
			("Centrality.ClusteringCoefficient",	"Node Centrality",	"Local Clustering Coefficient",
				True,	funcScores,	"Score",				centrality.LocalClusteringCoefficient,	(G, )),
			("Centrality.PageRank", 				"Node Centrality",	"PageRank",
				True,	funcScores,	"Score",				centrality.PageRank, 					(G, )),
			("Centrality.KPath", 					"Node Centrality",	"k-Path Centrality",
				True,	funcScores,	"Score",				centrality.KPathCentrality,				(G, )),
			("Centrality.Katz",						"Node Centrality",	"Katz Centrality",
				True,	funcScores,	"Score",				centrality.KatzCentrality,				(G, )),
			("Centrality.Betweenness", 				"Node Centrality",	"Betweenness",
				True,	funcScores,	"Score",				centrality.ApproxBetweenness2,			(G, max(42, math.log(G.numberOfNodes()), True))),
			("Centrality.Closeness",				"Node Centrality",	"Closeness",
				True,	funcScores,	"Score",				centrality.ApproxCloseness,				(G, max(42, math.log(G.numberOfNodes()), True))),
			("Partition.Communities", 				"Partition",		"Communities",
				False,	funcSizes,	"Nodes per Community",	community.PLM,			 				(G, )),
			("Partition.ConnectedComponents", 		"Partition",		"Connected Components",
				False,	funcSizes,	"Nodes per Component",	classConnectedComponents,				(G, )),
			("Partition.CoreDecomposition", 		"Partition",		"k-Core Decomposition",
				False,	funcSizes,	"Nodes per Shell",		centrality.CoreDecomposition, 			(G, ))
		]: result.__addMeasure(parameter)

		if cls.__verbose:
			timerAll = stopwatch.Timer()
		result.__loadProperties()
		result.__loadMeasures()
		if cls.__verbose:
			if cls.__verboseLevel < 1:
				print("")
			print("\ntotal time (measures + stats + correlations): {:.2F} s".format(timerAll.elapsed))
			print("total speed: {:.1F} edges/s".format(G.numberOfEdges() / timerAll.elapsed))
		return result;


	@classmethod
	def setVerbose(cls, verbose=False, level=0, filename=""):
		""" set verbose behaviour of all public methods 
		
		Args:
			verbose: enable/disable display verbose
			level: set level of verbose (0, 1)
			filename: enable/disable additional logfile support to given file
		"""
		cls.__verbose = verbose
		cls.__verboseLevel = level
		cls.__verboseFilename = filename


	@classmethod
	def getVerbose(cls):
		""" returns verbose settings """
		return (cls.__verbose, cls.__verboseLevel, cls.__verboseFilename)


	@classmethod
	def setParallel(cls, parallel):
		""" set the number of parallel threads/processes to use """
		if (parallel < 1):
			raise ValueError("parallel < 1");
		cls.__parallel = parallel


	@classmethod
	def getParallel(cls):
		""" returns the number of parallel threads/processes to use """
		return cls.__parallel


	def getStat(self, measure):
		""" returns an directory of all computed measure values """
		return self.__measures[measure]["stat"]


	def getCategory(self, measure):
		""" returns the category of an given measure """
		return self.__measures[measure]["category"]


	def getElapsedTime(self, measure):
		""" returns the elapsed computation times of an given measure """
		return self.__measures[measure]["time"]


	def output(self, outputType, directory, style="light", color=(0, 0, 1), parallel=False):
		""" outputs a computed profile to disk
		
		Args:
			outputType: profile output format ("HTML", "LaTeX")
			directory: directory to write
			style: style of generated output ("light")
			color: mainly used color of given style
			arallel: run some additional parts of the generation in parallel (experimental)
		"""
		
		# TODO: type -> enum
		options_type = ["HTML", "LaTeX", None]
		for o in options_type:
			if o is None:
				raise ValueError("unknown output type: options are " + str(options_type[0:len(options_type)-1]))
			if o == outputType:
				break

		filename  = "{0}.".format(self.__G.getName())

		result = self.__format(
			type = outputType,
			directory = directory,
			filename = filename,
			style = style,
			color = color,
			pageIndex = 0,
			parallel = parallel
		)

		if outputType == "HTML":
			js = readfile("html/profiling.js", False).replace("\\\\", "\\")
			css = readfile("html/profiling.css", False).replace("\\\\", "\\")
			result = """
				<!DOCTYPE HTML>
				<html>
					<head>
						<meta charset="utf-8">
						<style>
							""" + css + """
						</style>
						<script type="text/x-mathjax-config">
							MathJax.Hub.Config({tex2jax: {inlineMath: [['$','$']]}});
						</script>
						<script type="text/javascript" src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>
						<script>
							""" + js + """
						</script>
					</head>
					<body>
						""" + result + """
					</body>
				</html>
			"""
			filename += "html"
		elif outputType == "LaTeX":
			result = result
			if directory[-1] == "/":
				directory += filename[:-1]
			else:
				directory += "/" + filename[:-1]
			filename += "tex"
		else:
			raise Error("unknown output type")

		with open(directory + "/" + filename, 'w') as file:
			file.write(result)


	def show(self, style="light", color=(0, 0, 1), parallel=False):
		""" display computed profile
		
		Args:
			style: style of generated output ("light")
			color: mainly used color of given style
			parallel: run some additional parts of the generation in parallel (experimental)
		"""
		try:
			__IPYTHON__
		except:
			raise RuntimeError("this function cannot be used outside ipython notebook")

		result = self.__format(
			type = "HTML",
			directory = "",
			filename = "",
			style = style,
			color = color,
			pageIndex = self.__pageCount,
			parallel = parallel
		)

		display_html(HTML(result))
		self.__pageCount = self.__pageCount + 1


	def __format(self, type, directory, filename, style, color, pageIndex, parallel):
		""" layouts the profile	"""
		theme = plot.Theme()
		theme.set(style, color)

		# TODO: refactor and use decorator design pattern instead if an 3rd format is supported 
		if type == "HTML":
			plottype = "SVG"
			options = []
			templateProfile = readfile("html/profile.html", False)
			templateMeasure = readfile("html/measure.html", False)
		elif type == "LaTeX":
			plottype = "PDF"
			if directory[-1] == "/":
				output_dir = directory + filename[:-1]
			else:
				output_dir = directory + "/" + filename[:-1]
			if not os.path.isdir(output_dir):
				os.mkdir(output_dir)
			options = [output_dir, filename]
			templateProfile = readfile("latex/profile.tex", False)
			templateMeasure = readfile("latex/measure.tex", False)

		# generate plots
		pool = multiprocessing.ThreadPool(self.__parallel, parallel)
		for name in self.__measures:
			category = self.__measures[name]["category"]
			pool.put(
				plot.Measure(plottype, options, name, (
					0,
					self.__measures[name]["stat"],
					category,
					self.__measures[name]["name"],
					self.__measures[name]["label"],
					theme
				))
			)
			pool.put(
				plot.Measure(plottype, options, name, (
					1,
					self.__measures[name]["stat"],
					category,
					self.__measures[name]["name"],
					self.__measures[name]["label"],
					theme
				))
			)
			if category == "Partition":
				pool.put(
					plot.Measure(plottype, options, name, (
						2,
						self.__measures[name]["stat"],
						category,
						self.__measures[name]["name"],
						self.__measures[name]["label"],
						theme
					))
				)
		for category in self.__correlations:
			keyBList = []
			for keyA in self.__measures:
				if self.__measures[keyA]["category"] == category and self.__measures[keyA]["correlate"]:
					keyBList.append(keyA)
					for keyB in keyBList:
						if keyA != keyB:
							try:
								value = self.__correlations[category][keyA][keyB]
								keyFirst = keyA
								keySecond = keyB
							except:
								value = self.__correlations[category][keyB][keyA]
								keyFirst = keyB
								keySecond = keyA
							pool.put(
								plot.Scatter(plottype, options, keyFirst, (
									keySecond,
									self.__measures[keyFirst]["name"],
									self.__measures[keySecond]["name"],
									self.__measures[keyFirst]["label"],
									self.__measures[keySecond]["label"],
									self.__measures[keyFirst]["stat"],
									self.__measures[keySecond]["stat"],
									value["stat"],
									theme
								))
							)
		while pool.numberOfTasks() > 0:
			(taskType, name, data) = pool.get()
			try:
				category = self.__measures[name]["category"]

				if taskType == "Plot.Measure":
					(index, image) = data
					self.__verbosePrint("Plot.Measure", level=1)
					self.__measures[name]["image"][index] = image
				
				elif taskType == "Plot.Scatter":
					(nameB, image) = data
					self.__verbosePrint("Plot.Scatter (" + name + " <-> " + nameB + ")", level=1)
					self.__correlations[category][name][nameB]["image"] = image
			except Exception as e:
				self.__verbosePrint("Error (Post Processing): " + taskType + " - " + name, level=-1)
				self.__verbosePrint(str(e), level=-1)
		pool.join()

		# outline results
		results = {}
		for category in self.__correlations:
			results[category] = {}
			results[category]["Correlations"] = {}
			results[category]["Correlations"]["HeatMaps"] = ""
			results[category]["Correlations"]["ScatterPlots"] = ""
			results[category]["Measures"] = ""
			results[category]["Overview"] = ""

			def funcHeatMap(category, correlationName):
				result = ""
				keyBList = []

				if type == "HTML":
					result += "<div class=\"SubCategory HeatTable\" data-title=\"" + correlationName + "\">"
					for keyA in self.__measures:
						if self.__measures[keyA]["category"] == category and self.__measures[keyA]["correlate"]:
							keyBList.append(keyA)
							for keyB in keyBList:
								try:
									value = self.__correlations[category][keyA][keyB]
								except:
									value = self.__correlations[category][keyB][keyA]
								nameA = self.__measures[keyA]["name"]
								nameB = self.__measures[keyB]["name"]
								result += "<div class=\"HeatCell\" title=\"" + nameB + " - " + nameA + "\" data-heat=\"{:+.3F}\"></div>".format(value["stat"]["Value"][correlationName])
							result += "<div class=\"HeatCellName\">" + nameB + "</div><br>"
					result += "</div>"
				elif type == "LaTeX":
					i = 0
					n = len(self.__measures)
					result += "\\subsection{" + correlationName + "}\n"
					result += "\\begin{tabular}{" + ("l" * (n+1)) +"}"
					for keyA in self.__measures:
						if self.__measures[keyA]["category"] == category and self.__measures[keyA]["correlate"]:
							keyBList.append(keyA)
							for keyB in keyBList:
								try:
									value = self.__correlations[category][keyA][keyB]
								except:
									value = self.__correlations[category][keyB][keyA]
								nameA = self.__measures[keyA]["name"]
								nameB = self.__measures[keyB]["name"]
								result += "\\cellcolor{"
								result += "red" if value["stat"]["Value"][correlationName] > 0 else "blue"
								result += "!" + str(abs(value["stat"]["Value"][correlationName]) * 100)
								result += "}"
								result += "{:+.3F} & ".format(value["stat"]["Value"][correlationName])
							result += "\\multicolumn{" + str(n-i) + "}" + "{l}{" + nameB + "} \\\\"
						i += 1
					result += "\\end{tabular}"
				return result
			if self.__config.getMeasureCorrelation("Pearson"):
				results[category]["Correlations"]["HeatMaps"] += funcHeatMap(category, "Pearson's Correlation Coefficient")
			if self.__config.getMeasureCorrelation("Spearman"):
				results[category]["Correlations"]["HeatMaps"] += funcHeatMap(category, "Spearman's Rank Correlation Coefficient")
			if self.__config.getMeasureCorrelation("Fechner"):
				results[category]["Correlations"]["HeatMaps"] += funcHeatMap(category, "Fechner's Correlation Coefficient")

			def funcScatterPlot(category):
				result = ""

				keyBList = []
				for keyA in self.__measures:
					if self.__measures[keyA]["category"] == category and self.__measures[keyA]["correlate"]:
						keyBList.append(keyA)
						for keyB in keyBList:
							if keyA != keyB:
								try:
									value = self.__correlations[category][keyA][keyB]
								except:
									value = self.__correlations[category][keyB][keyA]
								if type == "HTML":
									result += "<div class=\"Thumbnail_ScatterPlot\" data-title=\"" + keyB + "\" data-title-second=\"" + keyA + "\"><img src=\"data:image/svg+xml;utf8," + value["image"] + "\" /></div>"
								elif type == "LaTeX":
									result += "\n\\includegraphics[width=0.5\\textwidth]{{\"" + value["image"] + "\"}.pdf}"
				return result
			results[category]["Correlations"]["ScatterPlots"] += funcScatterPlot(category)

		# group results to fit templates (measure)
		for key in self.__measures:
			measure = self.__measures[key]
			name = measure["name"]
			category = measure["category"]
			image = measure["image"]
			stat = measure["stat"]
			assortativity = measure["assortativity"]
			centralization = measure["centralization"]
			algorithm = measure["class"].__name__

			description = "N/A"
			try:
				description = readfile("description/" + key + ".txt")
			except:
				pass

			extentions = ""
			try:
				if type == "HTML":
					extentions = "<div class=\"PartitionPie\"><img src=\"data:image/svg+xml;utf8," + image[2] + "\" /></div>"
				elif type == "LaTeX":
					extentions = "\n\\includegraphics[width=0.5\\textwidth]{{\"" + image[2] + "\"}.pdf}"

			except:
				pass
				
			results[category]["Measures"] += self.__formatMeasureTemplate(
				templateMeasure,
				pageIndex,
				key,
				name,
				image,
				stat,
				assortativity,
				centralization,
				extentions,
				description,
				algorithm
			)
			if type == "HTML":
				results[category]["Overview"] += "<div class=\"Thumbnail_Overview\" data-title=\"" + name + "\"><a href=\"#NetworKit_Page_" + str(pageIndex) + "_" + key + "\"><img src=\"data:image/svg+xml;utf8," + image[1] + "\" /></a></div>"
			elif type == "LaTeX":
				results[category]["Overview"] += "\n\\includegraphics[width=0.5\\textwidth]{{\"" + image[1] + "\"}.pdf}\n"

		result = self.__formatProfileTemplate(
			templateProfile,
			pageIndex,
			results
		)
		return result


	def __formatMeasureTemplate(self, template, pageIndex, key, name, image, stat, assortativity, centralization, extentions, description, algorithm):
		""" format measure template - all function parameters, are available for the template """
		result = template.format(**locals())
		return result


	def __formatProfileTemplate(self, template, pageIndex, results):
		""" format profile template - all function parameters, are available for the template """
		properties = self.__properties
		result = template.format(**locals())
		return result


	def __addMeasure(self, args):
		""" add a measure if it is enabled in the given config object """
		(measureKey, measureCategory, measureName, correlate, getter, label, measureClass, parameters) = args
		if self.__config.getMeasure(measureKey):
			measure = {}
			measure["name"] = measureName
			measure["category"] = measureCategory
			measure["correlate"] = correlate
			measure["getter"] = getter
			measure["label"] = label
			measure["class"] = measureClass
			measure["parameters"] = parameters
			measure["data"] = {}
			measure["image"] = {}
			self.__measures[measureKey] = measure
		try:
			self.__correlations[measureCategory]
		except:
			self.__correlations[measureCategory] = {}


	def __loadProperties(self):
		""" calculate the network properties """
		self.__properties["Name"] = self.__G.getName()
		self.__properties["Nodes"] = self.__G.numberOfNodes()
		self.__properties["Edges"] = self.__G.numberOfEdges()
		self.__properties["Density"] = self.__G.density()
		self.__properties["Directed"] = self.__G.isDirected()
		self.__properties["Weighted"] = self.__G.isWeighted()
		self.__properties["Self Loops"] = self.__G.numberOfSelfLoops()

		if self.__config.getProperty("Diameter"):
			try:
				timerInstance = stopwatch.Timer()
				self.__verbosePrint("Diameter: ", end="")
				diameter = properties.Diameter.estimatedDiameterRange(self.__G, error=0.1)
				elapsedMain = timerInstance.elapsed
				self.__verbosePrint("{:.2F} s".format(elapsedMain))
				self.__verbosePrint("")
			except:
				self.__verbosePrint("Diameter raised exception")
				diameter = "N/A"
		else:
			diameter = "N/A"
		self.__properties["Diameter Range"] = diameter


		timerInstance = stopwatch.Timer()
		self.__verbosePrint("Connected Components: ", end="")
		try:
			if self.__G.isDirected():
				cc = components.StronglyConnectedComponents(self.__G)
			else:
				cc = components.ConnectedComponents(self.__G)
			cc.run()
			num = cc.numberOfComponents()
		except:
			self.__verbosePrint("ConnectedComponents raised exception")
			num = "N/A"
		elapsedMain = timerInstance.elapsed
		self.__verbosePrint("{:.2F} s".format(elapsedMain))
		self.__verbosePrint("")
		self.__properties["Connected Components"] = num


	def __loadMeasures(self):
		""" calculate the network measures and stats """
		pool = multiprocessing.ThreadPool(self.__parallel, False)

		for name in self.__measures:
			measure = self.__measures[name]
			self.__verbosePrint(name + ": ", end="")
			try:
				instance = measure["class"](*measure["parameters"])
			except Exception as e:
				del self.__measures[name]
				self.__verbosePrint("(removed)\n>> " + str(e))
				continue

			# run algorithm and get result
			timerInstance = stopwatch.Timer()
			instance.run()
			measure["data"]["sample"] = measure["getter"](instance)
			#self.__verbosePrint("{0} called on {1}".format(measure["getter"], instance))
			elapsedMain = timerInstance.elapsed
			self.__verbosePrint("{:.2F} s".format(elapsedMain))

			self.__verbosePrint("    Sort: ", end="")
			timerPostSort = stopwatch.Timer()
			measure["data"]["sorted"] = stat.sorted(measure["data"]["sample"])
			elapsedPostSort = timerPostSort.elapsed
			self.__verbosePrint("{:.2F} s".format(elapsedPostSort))

			self.__verbosePrint("    Rank: ", end="")
			timerPostRank = stopwatch.Timer()
			measure["data"]["ranked"] = stat.ranked(measure["data"]["sample"])
			elapsedPostRank = timerPostRank.elapsed
			self.__verbosePrint("{:.2F} s".format(elapsedPostRank))

			if self.__measures[name]["category"] == "Node Centrality":
				self.__verbosePrint("    Assortativity: ", end="")
				timerPostAssortativity = stopwatch.Timer()
				assortativity = properties.Assortativity(self.__G, measure["data"]["sample"])
				assortativity.run()
				measure["assortativity"] = assortativity.getCoefficient()
				elapsedPostAssortativity = timerPostAssortativity.elapsed
				self.__verbosePrint("{:.2F} s".format(elapsedPostAssortativity))
			else:
				measure["assortativity"] = float("nan")

			if self.__measures[name]["category"] == "Node Centrality":
				self.__verbosePrint("    Centralization: ", end="")
				timerPostCentralization = stopwatch.Timer()
				try:
					measure["centralization"] = instance.centralization()
				except:
					self.__verbosePrint("Centrality.centralization not properly defined for {0}. ".format(name), level=0, end="")
					measure["centralization"] = float("nan")
				elapsedPostCentralization = timerPostCentralization.elapsed
				self.__verbosePrint("{:.2F} s".format(elapsedPostCentralization))
			else:
				measure["centralization"] = float("nan")

			measure["time"] = (
				elapsedMain,
				elapsedPostSort,
				elapsedPostRank,
				elapsedPostAssortativity,
				elapsedPostCentralization
			)

		self.__verbosePrint("")

		for name in self.__measures:
			# the fix below avoids a cake-diagram for connected graphs,
			# as there is no other other way to include connected components for now,
			# the detailed output for connected graphs is accepted
			#if len(self.__measures[name]["data"]["sample"]) <= 1:	# for case connected graph...
			#	del self.__measures[name]
			#else:
			category = self.__measures[name]["category"]
			pool.put(
				stat.Stat(name, (
					self.__measures[name]["data"]["sample"],
					self.__measures[name]["data"]["sorted"],
					self.__measures[name]["data"]["ranked"],
					category == "Partition"
				))
			)

		while pool.numberOfTasks() > 0:
			(type, name, data) = pool.get()

			try:
				category = self.__measures[name]["category"]

				if type == "Stat":
					self.__measures[name]["stat"] = data
					self.__verbosePrint("Stat: " + name, level=1)
					if self.__measures[name]["correlate"]:
						for key in self.__correlations[category]:
							self.__correlations[category][key][name] = {}
							self.__correlations[category][key][name]["stat"] = {}
							pool.put(
								stat.Correlation(key, (
									name,
									self.__measures[key]["data"]["sample"],
									self.__measures[key]["data"]["ranked"],
									self.__measures[key]["stat"],
									self.__measures[name]["data"]["sample"],
									self.__measures[name]["data"]["ranked"],
									self.__measures[name]["stat"]
								))
							)
						self.__correlations[category][name] = {}
						self.__correlations[category][name][name] = {}
						self.__correlations[category][name][name]["stat"] = {}
						self.__correlations[category][name][name]["stat"]["Value"] = {
							"Spearman's Rank Correlation Coefficient": 1,
							"Pearson's Correlation Coefficient": 1,
							"Fechner's Correlation Coefficient": 1
						}
						self.__correlations[category][name][name]["image"] = ""

				elif type == "Correlation":
					(nameB, correlation) = data
					self.__verbosePrint("Correlation: " + name + " <-> " + nameB, level=1)
					self.__correlations[category][name][nameB]["stat"] = correlation
			except Exception as e:
				self.__verbosePrint("Error (Post Processing): " + type + " - " + nam-e, level=-1)
				self.__verbosePrint(str(e), level=-1)

		pool.join()


	@classmethod
	def __verbosePrint(cls, text="", end="\n", level=0):
		""" print for verbose output """
		if cls.__verboseLevel >= level:
			text = text + end
		else:
			text = "."

		if cls.__verbose or level < 0:
			print(text, end="", flush=True)

		if cls.__verboseFilename != "":
			with open(cls.__verboseFilename, 'a+') as file:
				file.write(text)