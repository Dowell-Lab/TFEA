from bokeh.charts import Bar,BoxPlot, HeatMap, output_file, show
from bokeh.models import HoverTool
import bokeh.layouts as layouts
import matplotlib.pyplot as plt
from matplotlib import cm
from bokeh.sampledata.autompg import autompg as df
import math
import sys

def bar_color(bardict,output=True):
	if output:
		# output to static HTML file
		output_file("results.html")

	# Create bar plot 
	hover = HoverTool(tooltips=[("gene","@genes")])
	bar = Bar(bardict,values='values', label='genes', title="Ranked Gene Biomarkers", xlabel = "Genes", ylabel = "Importance Score", color='pval',legend=False,tools=[hover])

	if output:
		show(bar)	

	return bar

def bar_nocolor(bardict,title,output=True):
	if output:
		output_file("results.html")
	bar = Bar(bardict,values='values',label='genes',title=title,xlabel="",ylabel="",legend=False)

def profile(boxplotdict,output=True):
	if output:
		output_file("results.html")
	
	boxplot = BoxPlot(boxplotdict, values='zscore',label='gene',color='sensitivity',group='sensitivity',title = "Cell Profile")
	
	return boxplot

def heatmap(heatmapdict,output=True):
	if output:
		output_file("results.html")

	heatmap = HeatMap(heatmapdict, y='year', x='fruit', values='fruit_count', stat=None, legend=False)
	
	if output:
		show(heatmap)	

	return heatmap

if __name__ == "__main__":
	bardict = {'values' : [1,2,3,4,5], 'genes' : ['a', 'b', 'c', 'd', 'e'], 'pval' : [.1,.2,.3,.4,.5]}
	output_file("results.html")
	bar = bar_color(bardict)
	heatmapdict = {'fruit': ['apples', 'apples', 'apples', 'apples', 'apples',
                    'pears', 'pears', 'pears', 'pears', 'pears',
                    'bananas', 'bananas', 'bananas', 'bananas', 'bananas'],
          'fruit_count': [4, 5, 8, 12, 4, 6, 5, 4, 8, 7, 1, 2, 4, 8, 12],
          'year': [2009, 2010, 2011, 2012, 2013, 2009, 2010, 2011, 2012, 2013, 2009, 2010,
                   2011, 2012, 2013]}

	boxplotdict = {'sensitivity' : ["non-sensitive","sensitive","non-sensitive","sensitive","non-sensitive","sensitive","non-sensitive","sensitive","non-sensitive","sensitive","non-sensitive","sensitive"], 'gene' : ['a','a','a','a','a','a','b','b','b','b','b','b'], 'zscore' : [1,10,0,12,2,13,1,20,3,40,5,40]}
	heatmap = heatmap(heatmapdict)
	boxplot = profile(boxplotdict)
	show(layouts.row(bar,boxplot))
