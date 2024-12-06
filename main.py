from fastapi import FastAPI, Request
from fastapi.templating import Jinja2Templates
from fastapi.staticfiles import StaticFiles
from fastapi.responses import HTMLResponse
from pydantic import BaseModel
from api.engine.diagnosis_assistant import DiagnosisAssistantEngine
import uvicorn
import os
from dotenv import load_dotenv
import sys

def verify_environment():
    """Verify that all required environment variables are set and valid"""
    load_dotenv()
    
    api_key = os.getenv('OPENAI_API_KEY')
    
    if not api_key:
        print("\033[91mError: OPENAI_API_KEY not found in .env file\033[0m")
        print("\033[93m1. Create a .env file in the root directory")
        print("2. Add your OpenAI API key like this: OPENAI_API_KEY=your-key-here")
        print("3. Make sure your API key is valid\033[0m")
        return False
        
    if not api_key.startswith('sk-') or len(api_key) < 20:
        print("\033[91mError: OPENAI_API_KEY appears to be invalid\033[0m")
        print("\033[93mMake sure you've entered a valid OpenAI API key in the .env file\033[0m")
        return False
        
    return True

app = FastAPI()

# Mount static files
app.mount("/static", StaticFiles(directory="templates/static"), name="static")

# Templates
templates = Jinja2Templates(directory="templates")

class DiagnosisRequest(BaseModel):
    symptoms: str

class DiagnosisResponse(BaseModel):
    diagnosis: str

engine = DiagnosisAssistantEngine()

@app.get("/", response_class=HTMLResponse)
async def home(request: Request):
    return templates.TemplateResponse("index.html", {"request": request})

@app.post("/diagnosis")
async def get_diagnosis(request: DiagnosisRequest) -> DiagnosisResponse:
    diagnosis = engine.generate_diagnosis(symptoms=request.symptoms)
    return DiagnosisResponse(diagnosis=diagnosis)

if __name__ == "__main__":
    if not verify_environment():
        sys.exit(1)
    
    print("\033[92m✓ Environment validated successfully\033[0m")
    print("\033[94mStarting server...\033[0m")
    uvicorn.run(app, host="0.0.0.0", port=8000)
