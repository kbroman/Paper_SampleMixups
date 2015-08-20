html: sample_mixups.html

sample_mixups.html: sample_mixups.asciidoc
	echo 'library(ascii); Asciidoc("sample_mixups.asciidoc", width=9, height=4, prefix.string="Figs/fig", driver=weaverAsciidoc(), pdf=FALSE)' | R --no-save --no-restore
	asciidoc -a data-uri -a toc -a toclevels=4 -a num sample_mixups.asciidoc.txt
	mv sample_mixups.asciidoc.html sample_mixups.html

adhara:	sample_mixups.html
	scp sample_mixups.html sample_mixups.asciidoc adhara.biostat.wisc.edu:public_html/attiedata/
	cp *.html ~/DropBox/eQTL\ Paper/

all: adhara

clean: 
	\rm -r *~ *.bak *.dvi *.log *.aux *.tex *.html *.asciidoc.txt
