# PySIO AI - API de Diagnóstico Médico

Una API inteligente para diagnóstico médico basado en síntomas utilizando inteligencia artificial y evidencia médica científica.

## 🚀 Características

- **Diagnóstico Inteligente**: Genera diagnósticos médicos basados en síntomas usando GPT-4
- **Evidencia Científica**: Integra búsquedas en PubMed para respaldar diagnósticos
- **Base de Conocimiento Vectorial**: Almacena y recupera información médica relevante
- **API RESTful**: Endpoints bien documentados y fáciles de usar
- **Arquitectura Limpia**: Separación clara entre modelos, servicios y controladores

## 🏗️ Arquitectura

```
app/
├── config/          # Configuración y variables de entorno
├── models/          # Modelos de datos Pydantic
├── services/        # Lógica de negocio y servicios
├── controllers/     # Endpoints de la API
└── utils/           # Utilidades y helpers
```

## 📋 Prerrequisitos

- Python 3.13+
- OpenAI API Key
- Email válido para PubMed/Entrez

## ⚡ Instalación Rápida

1. **Clonar el repositorio**

   ```bash
   git clone <repository-url>
   cd pysio_ai
   ```

2. **Crear entorno virtual**

   ```bash
   python3 -m venv env
   source env/bin/activate  # En macOS/Linux
   # o
   env\Scripts\activate     # En Windows
   ```

3. **Instalar dependencias**

   ```bash
   pip install -r requirements.txt
   ```

4. **Configurar variables de entorno**

   ```bash
   cp .env.example .env
   # Editar .env con tu API key de OpenAI
   ```

5. **Ejecutar la API**
   ```bash
   python run.py
   ```

## 🔧 Uso del Makefile

```bash
# Ver todos los comandos disponibles
make help

# Instalar dependencias
make install

# Ejecutar en modo desarrollo
make run-dev

# Ejecutar en modo producción
make run-prod

# Ejecutar tests
make test

# Verificar configuración
make check-env
```

## 📚 Endpoints de la API

### Diagnóstico

- `POST /api/v1/diagnosis` - Generar diagnóstico médico
- `GET /api/v1/diagnosis/history` - Obtener historial de diagnósticos
- `POST /api/v1/diagnosis/update-knowledge` - Actualizar base de conocimiento

### Utilidades

- `GET /api/v1/health` - Verificar estado de la API
- `GET /docs` - Documentación interactiva (Swagger UI)
- `GET /redoc` - Documentación alternativa

## 🔍 Ejemplo de Uso

### Generar un Diagnóstico

```bash
curl -X POST "http://localhost:8000/api/v1/diagnosis" \
     -H "Content-Type: application/json" \
     -d '{
       "symptoms": [
         {
           "description": "Dolor de cabeza intenso",
           "severity": "moderate",
           "duration": "2 hours"
         }
       ],
       "patient_age": 35,
       "patient_gender": "female",
       "medical_history": "Sin antecedentes relevantes"
     }'
```

### Respuesta

```json
{
  "diagnosis": "Migraña tensional",
  "confidence": 0.85,
  "recommendations": [
    "Descansar en un ambiente tranquilo",
    "Aplicar compresas frías",
    "Considerar analgésicos de venta libre"
  ],
  "timestamp": "2024-01-15T10:30:00Z"
}
```

## ⚙️ Configuración

### Variables de Entorno

| Variable             | Descripción                 | Valor por Defecto |
| -------------------- | --------------------------- | ----------------- |
| `OPENAI_API_KEY`     | API Key de OpenAI           | Requerido         |
| `OPENAI_MODEL`       | Modelo de OpenAI a usar     | `gpt-4o-mini`     |
| `OPENAI_TEMPERATURE` | Temperatura para generación | `0.7`             |
| `HOST`               | Host del servidor           | `0.0.0.0`         |
| `PORT`               | Puerto del servidor         | `8000`            |
| `DEBUG`              | Modo debug                  | `false`           |
| `ENTREZ_EMAIL`       | Email para PubMed           | Requerido         |
| `LOG_LEVEL`          | Nivel de logging            | `INFO`            |

## 🧪 Testing

```bash
# Ejecutar todos los tests
make test

# Ejecutar tests con coverage
make test

# Ejecutar tests en modo watch
make test-watch
```

## 🚀 Despliegue

### Desarrollo

```bash
make run-dev
```

### Producción

```bash
make run-prod
```

### Docker

```bash
make docker-build
make docker-run
```

## 📊 Monitoreo

- **Health Check**: `/api/v1/health`
- **Logs**: Estructurados en formato JSON
- **Métricas**: Endpoints de estado y rendimiento

## 🤝 Contribución

1. Fork el proyecto
2. Crea una rama para tu feature (`git checkout -b feature/AmazingFeature`)
3. Commit tus cambios (`git commit -m 'Add some AmazingFeature'`)
4. Push a la rama (`git push origin feature/AmazingFeature`)
5. Abre un Pull Request

## 📝 Licencia

Este proyecto está bajo la Licencia MIT. Ver el archivo `LICENSE` para más detalles.

## ⚠️ Descargo de Responsabilidad

**IMPORTANTE**: Esta API es solo para fines educativos y de investigación. No debe usarse para diagnóstico médico real. Siempre consulta con profesionales de la salud calificados para cualquier problema médico.

## 🆘 Soporte

- **Issues**: [GitHub Issues](https://github.com/your-repo/issues)
- **Documentación**: `/docs` endpoint en la API
- **Email**: [tu-email@example.com]

---

Desarrollado con ❤️ para la comunidad médica y de IA
