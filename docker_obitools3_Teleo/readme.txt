
Hier wordt een docker image gemaakt om obitools3 te runnen. Het welkomstscript wordt eerst gerund en dan wordt de container interactief geopend.

In docker voor windows geïnstalleerd is, kan je gewoon een windows powershell openen in deze directory.

#Build image (naam obitools3pv)
docker build -t obitools3pv .

#Run image (Het Dockerfile bestand regelt dat de image in interactieve modus open is)
docker run --rm -it obitools3pv

#Run met gekende naam (Let op, naam moet uniek zijn)
docker run -it --name obitools3testcontainer obitools3pv

#Script uitvoeren in container (powershell)
docker cp test.sh obitools3testcontainer:/app/test.sh
docker stop obitools3testcontainer
docker start obitools3testcontainer