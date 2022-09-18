GMSH:=gmsh
FORMAT?=msh
SIZE?=0.01
DATA=data
CASE?=1
NAME=$(DATA)/mesh$(CASE).$(FORMAT)

all:
	@if [ ! -d $(DATA) ] ; \
	then \
		mkdir $(DATA) ; \
	fi;
	$(GMSH) mesh$(CASE).geo -2 -format $(FORMAT) -setnumber h $(SIZE) -o $(NAME)
	python3 mshtoxdmf$(CASE).py

clean:
	rm -f -r data/
