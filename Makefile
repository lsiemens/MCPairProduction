build:
	sphinx-apidoc -f -o ./doc/source/apidoc ./

deploy:
	cp -r ./doc/build/html/* ../MCPairProduction-pages/
