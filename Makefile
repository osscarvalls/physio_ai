# Makefile para PySIO AI - Aplicación Refactorizada
# Arquitectura: Models, Services, Controllers

.PHONY: help install run test clean lint format docs

# Variables
PYTHON = python3
VENV = env
PIP = $(VENV)/bin/pip
PYTHON_VENV = $(VENV)/bin/python
APP_NAME = pysio_ai

# Colores para output
GREEN = \033[92m
BLUE = \033[94m
RED = \033[91m
YELLOW = \033[93m
RESET = \033[0m

help: ## Muestra esta ayuda
	@echo "$(BLUE)PySIO AI - Makefile$(RESET)"
	@echo "$(GREEN)Comandos disponibles:$(RESET)"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  $(YELLOW)%-15s$(RESET) %s\n", $$1, $$2}'

install: ## Instala las dependencias del proyecto
	@echo "$(BLUE)Instalando dependencias...$(RESET)"
	$(PYTHON) -m venv $(VENV)
	$(PIP) install --upgrade pip
	$(PIP) install -r requirements.txt
	@echo "$(GREEN)✅ Dependencias instaladas correctamente$(RESET)"

install-dev: ## Instala dependencias de desarrollo
	@echo "$(BLUE)Instalando dependencias de desarrollo...$(RESET)"
	$(PIP) install -r requirements.txt
	$(PIP) install pytest pytest-cov pytest-asyncio black isort mypy
	@echo "$(GREEN)✅ Dependencias de desarrollo instaladas$(RESET)"

run: ## Ejecuta la aplicación
	@echo "$(BLUE)🚀 Iniciando PySIO AI...$(RESET)"
	$(PYTHON_VENV) run.py

run-dev: ## Ejecuta la aplicación en modo desarrollo
	@echo "$(BLUE)🔧 Iniciando PySIO AI en modo desarrollo...$(RESET)"
	$(PYTHON_VENV) -m uvicorn app.main:app --host 0.0.0.0 --port 8000 --reload

run-prod: ## Ejecuta la aplicación en modo producción
	@echo "$(BLUE)🏭 Iniciando PySIO AI en modo producción...$(RESET)"
	$(PYTHON_VENV) -m uvicorn app.main:app --host 0.0.0.0 --port 8000 --workers 4

test: ## Ejecuta los tests
	@echo "$(BLUE)🧪 Ejecutando tests...$(RESET)"
	$(PYTHON_VENV) -m pytest app/tests/ -v --cov=app --cov-report=html

test-watch: ## Ejecuta tests en modo watch
	@echo "$(BLUE)👀 Ejecutando tests en modo watch...$(RESET)"
	$(PYTHON_VENV) -m pytest app/tests/ -v --cov=app -f

lint: ## Ejecuta el linter
	@echo "$(BLUE)🔍 Ejecutando linter...$(RESET)"
	$(PYTHON_VENV) -m flake8 app/ --max-line-length=100 --exclude=__pycache__

format: ## Formatea el código
	@echo "$(BLUE)🎨 Formateando código...$(RESET)"
	$(PYTHON_VENV) -m black app/ --line-length=100
	$(PYTHON_VENV) -m isort app/ --profile=black

format-check: ## Verifica el formato del código
	@echo "$(BLUE)🔍 Verificando formato del código...$(RESET)"
	$(PYTHON_VENV) -m black app/ --check --line-length=100
	$(PYTHON_VENV) -m isort app/ --check-only --profile=black

type-check: ## Verifica tipos con mypy
	@echo "$(BLUE)🔍 Verificando tipos...$(RESET)"
	$(PYTHON_VENV) -m mypy app/ --ignore-missing-imports

clean: ## Limpia archivos temporales
	@echo "$(BLUE)🧹 Limpiando archivos temporales...$(RESET)"
	find . -type f -name "*.pyc" -delete
	find . -type d -name "__pycache__" -delete
	find . -type d -name "*.egg-info" -delete
	find . -type d -name ".pytest_cache" -delete
	find . -type d -name "htmlcov" -delete
	rm -rf .coverage
	rm -rf build/
	rm -rf dist/
	@echo "$(GREEN)✅ Limpieza completada$(RESET)"

clean-venv: ## Elimina el entorno virtual
	@echo "$(RED)🗑️ Eliminando entorno virtual...$(RESET)"
	rm -rf $(VENV)
	@echo "$(GREEN)✅ Entorno virtual eliminado$(RESET)"

docs: ## Genera documentación
	@echo "$(BLUE)📚 Generando documentación...$(RESET)"
	$(PYTHON_VENV) -m pdoc app/ --html --output-dir docs/
	@echo "$(GREEN)✅ Documentación generada en docs/$(RESET)"

setup: ## Configuración inicial del proyecto
	@echo "$(BLUE)⚙️ Configurando proyecto...$(RESET)"
	@if [ ! -f .env ]; then \
		echo "$(YELLOW)Creando archivo .env...$(RESET)"; \
		cp .env.example .env 2>/dev/null || echo "OPENAI_API_KEY=your-key-here" > .env; \
		echo "$(GREEN)✅ Archivo .env creado. Por favor, configura tu API key.$(RESET)"; \
	else \
		echo "$(GREEN)✅ Archivo .env ya existe$(RESET)"; \
	fi
	@if [ ! -d docs ]; then mkdir -p docs; fi
	@if [ ! -d app/tests ]; then mkdir -p app/tests; fi
	@echo "$(GREEN)✅ Configuración completada$(RESET)"

check-env: ## Verifica la configuración del entorno
	@echo "$(BLUE)🔍 Verificando configuración del entorno...$(RESET)"
	$(PYTHON_VENV) -c "from app.config.settings import settings; print('✅ Configuración válida')"
	@echo "$(GREEN)✅ Entorno configurado correctamente$(RESET)"

docker-build: ## Construye imagen Docker
	@echo "$(BLUE)🐳 Construyendo imagen Docker...$(RESET)"
	docker build -t $(APP_NAME) .

docker-run: ## Ejecuta contenedor Docker
	@echo "$(BLUE)🚀 Ejecutando contenedor Docker...$(RESET)"
	docker run -p 8000:8000 --env-file .env $(APP_NAME)

# Comando por defecto
.DEFAULT_GOAL := help 