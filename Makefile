# Variables
PYTHON_VERSION = 3.11.5
VENV_NAME = env
PYTHON = $(VENV_NAME)/bin/python
PIP = $(VENV_NAME)/bin/pip

# Colors for terminal output
CYAN = \033[0;36m
GREEN = \033[0;32m
YELLOW = \033[0;33m
RED = \033[0;31m
NC = \033[0m # No Color

.PHONY: all install clean setup-pyenv setup-venv install-deps help

# Default target
all: help

help:
	@echo "$(CYAN)Available commands:$(NC)"
	@echo "$(GREEN)make install$(NC)    - Set up everything (pyenv, virtualenv, and dependencies)"
	@echo "$(GREEN)make clean$(NC)      - Remove virtual environment and cached files"
	@echo "$(GREEN)make setup-pyenv$(NC) - Install Python $(PYTHON_VERSION) using pyenv"
	@echo "$(GREEN)make setup-venv$(NC)  - Create and activate virtual environment"
	@echo "$(GREEN)make install-deps$(NC)- Install project dependencies"

install: setup-pyenv setup-venv install-deps
	@echo "$(GREEN)Installation complete!$(NC)"
	@echo "$(YELLOW)To activate the virtual environment, run:$(NC)"
	@echo "source $(VENV_NAME)/bin/activate"

setup-pyenv:
	@echo "$(CYAN)Setting up Python $(PYTHON_VERSION) with pyenv...$(NC)"
	@if ! command -v pyenv >/dev/null 2>&1; then \
		echo "$(YELLOW)Installing pyenv...$(NC)"; \
		brew install pyenv; \
	fi
	@if ! pyenv versions | grep $(PYTHON_VERSION) >/dev/null 2>&1; then \
		echo "$(YELLOW)Installing Python $(PYTHON_VERSION)...$(NC)"; \
		pyenv install $(PYTHON_VERSION); \
	fi
	@pyenv local $(PYTHON_VERSION)
	@echo "$(GREEN)Python $(PYTHON_VERSION) is set up!$(NC)"

setup-venv:
	@echo "$(CYAN)Setting up virtual environment...$(NC)"
	@if [ ! -d "$(VENV_NAME)" ]; then \
		python -m venv $(VENV_NAME); \
		echo "$(GREEN)Virtual environment created!$(NC)"; \
	else \
		echo "$(YELLOW)Virtual environment already exists.$(NC)"; \
	fi
	@$(PIP) install --upgrade pip

install-deps:
	@echo "$(CYAN)Installing dependencies...$(NC)"
	@$(PIP) install fastapi uvicorn python-dotenv jinja2 python-multipart openai
	@if [ -f "requirements.txt" ]; then \
		$(PIP) install -r requirements.txt; \
	fi

clean:
	@echo "$(CYAN)Cleaning up...$(NC)"
	@rm -rf $(VENV_NAME)
	@find . -type d -name "__pycache__" -exec rm -rf {} +
	@find . -type f -name "*.pyc" -delete
	@echo "$(GREEN)Cleanup complete!$(NC)" 