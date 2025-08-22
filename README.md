# 🏥 PySIO AI - Asistente de Diagnóstico Médico

[![Python](https://img.shields.io/badge/Python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![FastAPI](https://img.shields.io/badge/FastAPI-0.104+-green.svg)](https://fastapi.tiangolo.com/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**PySIO AI** es un sistema inteligente de asistencia médica que utiliza **Inteligencia Artificial** para generar diagnósticos preliminares basados en síntomas del paciente y evidencia médica de **PubMed**.

## 🚀 Características Principales

- **🤖 IA Avanzada**: Utiliza GPT-4 para análisis de síntomas
- **📚 Base de Conocimiento**: Integración con PubMed para evidencia médica
- **🏗️ Arquitectura Limpia**: Separación clara en Models, Services y Controllers
- **🔍 Búsqueda Semántica**: Recuperación inteligente de información médica
- **📊 API REST**: Endpoints bien documentados con FastAPI
- **⚡ Asíncrono**: Operaciones no bloqueantes para mejor rendimiento
- **🔒 Seguro**: Validación de datos y manejo de errores robusto

## 🏛️ Arquitectura del Sistema

```
┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
│   Controllers   │    │     Services    │    │      Models     │
│                 │    │                 │    │                 │
│ • API Endpoints │───▶│ • Lógica de     │───▶│ • Estructura    │
│ • Validación    │    │   Negocio       │    │   de Datos      │
│ • Manejo HTTP   │    │ • Integración   │    │ • Validación    │
└─────────────────┘    │   Externa       │    │ • Serialización │
                       └─────────────────┘    └─────────────────┘
```

### **Capas de la Arquitectura**

- **🎮 Controllers**: Manejan requests HTTP y respuestas de la API
- **🔧 Services**: Implementan la lógica de negocio y coordinación
- **📊 Models**: Definen la estructura de datos con validación Pydantic

## 📦 Instalación

### **Requisitos Previos**

- Python 3.11+
- pip
- make (opcional, para usar el Makefile)

### **Instalación Rápida**

```bash
# 1. Clonar el repositorio
git clone https://github.com/tu-usuario/pysio_ai.git
cd pysio_ai

# 2. Crear entorno virtual e instalar dependencias
make install

# 3. Configurar variables de entorno
cp .env.example .env
# Editar .env con tu API key de OpenAI

# 4. Verificar instalación
make check-env

# 5. Ejecutar la aplicación
make run
```

### **Instalación Manual**

```bash
# Crear entorno virtual
python3 -m venv venv
source venv/bin/activate  # En Windows: venv\Scripts\activate

# Instalar dependencias
pip install -r requirements.txt

# Configurar .env
cp .env.example .env
# Editar .env con tu configuración

# Ejecutar
python run.py
```

## ⚙️ Configuración

### **Variables de Entorno**

Crea un archivo `.env` en el directorio raíz:

```bash
# OpenAI Configuration
OPENAI_API_KEY=your-openai-api-key-here
OPENAI_ORGANIZATION_ID=your-org-id-here
OPENAI_MODEL=gpt-4o-mini
OPENAI_TEMPERATURE=0.0

# Server Configuration
HOST=0.0.0.0
PORT=8000
DEBUG=false

# PubMed Configuration
ENTREZ_EMAIL=your-email@example.com

# File Paths
STATIC_DIR=templates/static
TEMPLATES_DIR=templates
DOCS_DIR=docs
CHROMA_PERSIST_DIR=./chroma_db
```

## 🚀 Uso

### **Iniciar la Aplicación**

```bash
# Modo desarrollo (con auto-reload)
make run-dev

# Modo producción
make run-prod

# Script directo
python run.py
```

### **Acceder a la API**

- **🌐 Aplicación Web**: http://localhost:8000
- **📚 Documentación API**: http://localhost:8000/docs
- **📖 ReDoc**: http://localhost:8000/redoc
- **❤️ Health Check**: http://localhost:8000/api/v1/health

### **Endpoints Principales**

#### **Diagnóstico Médico**

```http
POST /api/v1/diagnosis
Content-Type: application/json

{
  "symptoms": "Dolor de cabeza intenso, náuseas, sensibilidad a la luz",
  "patient_age": 35,
  "patient_gender": "femenino",
  "medical_history": "Migrañas ocasionales"
}
```

#### **Respuesta**

```json
{
  "diagnosis": "Basado en los síntomas descritos, se sugiere una migraña...",
  "confidence": 0.85,
  "recommendations": [
    "Consulte con un profesional médico para confirmar el diagnóstico",
    "Mantenga un registro de sus síntomas y su evolución"
  ],
  "timestamp": "2024-01-15T10:30:00"
}
```

## 🛠️ Comandos del Makefile

```bash
# Ayuda
make help

# Instalación
make install          # Instala dependencias básicas
make install-dev      # Instala dependencias de desarrollo

# Ejecución
make run             # Ejecuta la aplicación
make run-dev         # Modo desarrollo con auto-reload
make run-prod        # Modo producción

# Testing y Calidad
make test            # Ejecuta tests
make test-watch      # Tests en modo watch
make lint            # Ejecuta linter
make format          # Formatea código
make type-check      # Verifica tipos

# Utilidades
make clean           # Limpia archivos temporales
make docs            # Genera documentación
make check-env       # Verifica configuración
```

## 🧪 Testing

```bash
# Ejecutar todos los tests
make test

# Tests con coverage
make test

# Tests en modo watch
make test-watch
```

## 📁 Estructura del Proyecto

```
pysio_ai/
├── app/                          # 🎯 Paquete principal
│   ├── __init__.py
│   ├── main.py                   # 🚀 Punto de entrada
│   ├── config/                   # ⚙️ Configuración
│   │   ├── __init__.py
│   │   └── settings.py           # Configuración centralizada
│   ├── models/                   # 📊 Modelos de datos
│   │   ├── __init__.py
│   │   └── diagnosis.py          # Modelos Pydantic
│   ├── services/                 # 🔧 Lógica de negocio
│   │   ├── __init__.py
│   │   ├── diagnosis_service.py  # Servicio principal
│   │   ├── llm_service.py        # Servicio LLM
│   │   └── evidence_service.py   # Servicio evidencia
│   ├── controllers/              # 🎮 Controladores API
│   │   ├── __init__.py
│   │   └── diagnosis_controller.py
│   └── utils/                    # 🛠️ Utilidades
├── templates/                     # 🎨 Templates HTML
├── docs/                         # 📚 Documentación
├── run.py                        # 🚀 Script de inicio
├── requirements.txt              # 📦 Dependencias
├── Makefile                      # 🛠️ Comandos de desarrollo
└── README.md                     # 📖 Este archivo
```

## 🔧 Desarrollo

### **Agregar Nuevas Funcionalidades**

1. **Modelos**: Definir en `app/models/`
2. **Servicios**: Implementar en `app/services/`
3. **Controladores**: Crear en `app/controllers/`
4. **Tests**: Agregar en `app/tests/`

### **Convenciones de Código**

- **Type Hints**: Usar en todas las funciones
- **Docstrings**: Documentar todas las clases y métodos
- **Logging**: Usar logging estructurado
- **Error Handling**: Manejar errores apropiadamente
- **Async/Await**: Usar para operaciones I/O

## 🚧 Roadmap

### **Versión 1.1** 🎯

- [ ] Tests unitarios completos
- [ ] Logging estructurado con structlog
- [ ] Cache para respuestas frecuentes
- [ ] Métricas y monitoreo

### **Versión 1.2** 🚀

- [ ] Autenticación y autorización
- [ ] Base de datos para historial
- [ ] API para múltiples idiomas
- [ ] Dashboard de administración

### **Versión 2.0** 🌟

- [ ] Machine Learning personalizado
- [ ] Integración con sistemas médicos
- [ ] Análisis de imágenes médicas
- [ ] Predicción de enfermedades

## 🤝 Contribución

¡Las contribuciones son bienvenidas! Por favor:

1. **Fork** el proyecto
2. **Crea** una rama para tu feature (`git checkout -b feature/AmazingFeature`)
3. **Commit** tus cambios (`git commit -m 'Add some AmazingFeature'`)
4. **Push** a la rama (`git push origin feature/AmazingFeature`)
5. **Abre** un Pull Request

### **Directrices de Contribución**

- Mantener la separación de responsabilidades
- Agregar tests para nuevas funcionalidades
- Documentar cambios en la API
- Seguir las convenciones de nomenclatura
- Usar type hints en todas las funciones

## 📄 Licencia

Este proyecto está bajo la Licencia MIT. Ver el archivo [LICENSE](LICENSE) para más detalles.

## ⚠️ Descargo de Responsabilidad

**IMPORTANTE**: PySIO AI es un **asistente médico** y **NO reemplaza** la consulta con profesionales médicos calificados. Siempre consulte con un médico para diagnósticos y tratamientos.

## 📞 Contacto

- **Proyecto**: [GitHub Issues](https://github.com/tu-usuario/pysio_ai/issues)
- **Email**: tu-email@example.com
- **Documentación**: [Wiki del Proyecto](https://github.com/tu-usuario/pysio_ai/wiki)

---

**🏥 PySIO AI** - Transformando el diagnóstico médico con Inteligencia Artificial

_Construido con ❤️ para la comunidad médica_
