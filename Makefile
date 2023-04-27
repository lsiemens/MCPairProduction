# Flags for apidoc
FLAGS = --separate -d=1 --no-toc --module-first

apidoc:
	sphinx-apidoc $(FLAGS) -f -o ./doc/source/apidoc ./

html:
	cd ./doc && make html

clean:
	cd ./doc && make clean

view:
	xdg-open ./doc/build/html/index.html

deploy:
	cp -r ./doc/build/html/* ../MCPairProduction-pages/
