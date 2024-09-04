from api.engine.diagnosis_assistant import DiagnosisAssistantEngine
from fastapi import FastAPI
from pydantic import BaseModel
import uvicorn


app = FastAPI()

class DiagnosisRequest(BaseModel):
    symptoms: str

class DiagnosisResponse(BaseModel):
    diagnosis: str

engine = DiagnosisAssistantEngine()

@app.post("/diagnosis")
async def get_diagnosis(request: DiagnosisRequest) -> DiagnosisResponse:
    # This is a placeholder implementation
    # In a real-world scenario, you would integrate with a medical diagnosis system
    diagnosis = engine.generate_diagnosis(symptoms=request)
    return DiagnosisResponse(diagnosis=diagnosis)

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000)
