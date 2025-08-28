# Makefile simple para PySIO AI
.PHONY: build run stop clean

build:
	sudo docker pull python:3.12-slim
	sudo docker build -t physio-ai .

run:
	sudo docker run -d --name physio-ai -p 8000:8000 --env-file .env physio-ai

stop:
	sudo docker stop physio-ai && sudo docker rm physio-ai

clean:
	sudo docker stop physio-ai 2>/dev/null || true
	sudo docker rm physio-ai 2>/dev/null || true


