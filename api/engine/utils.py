import os

from dotenv import load_dotenv


def get_prompts():
        prompts = {}
        prompts_dir = './api/prompts'
        for filename in os.listdir(prompts_dir):
            if filename.endswith('.txt'):
                file_path = os.path.join(prompts_dir, filename)
                with open(file_path, 'r') as file:
                    content = file.read()
                    prompts[filename[:-4]] = content  # Remove '.txt' from the key
        return prompts

def set_openai():
    load_dotenv()
    OPENAI_API_KEY = os.getenv('OPENAI_API_KEY')
    os.environ['OPENAI_API_KEY'] = OPENAI_API_KEY