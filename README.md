# PhysioAI - Physiotherapy Diagnosis Assistant

PhysioAI is an AI-powered physiotherapy diagnosis assistant that helps analyze symptoms and provide evidence-based medical insights. It combines natural language processing with medical knowledge to assist in preliminary diagnosis and clinic assessment.

## 🚨 Important Disclaimer

This tool is for educational and research purposes only. It should not be used as a substitute for professional medical advice, diagnosis, or treatment. Always seek the advice of your physician or other qualified health provider with any questions you may have regarding a medical condition.

## 🌟 Features

- **Symptom Analysis**: Process natural language descriptions of symptoms using OpenAI GPT models and LangChain
- **Evidence-Based Insights**: Retrieves relevant medical literature from PubMed using BioPython and XMLtodict
- **Modern Web Interface**: Clean and responsive design
- **Real-time Processing**: Immediate feedback and analysis through asynchronous FastAPI endpoints
- **Secure**: Requires OpenAI API key authentication managed via python-dotenv

## 🛠 Prerequisites

Before you begin, ensure you have:
- macOS or Linux operating system
- Python 3.11+ installed (via pyenv or system package manager)
- OpenAI API key (get one at https://platform.openai.com)

## 🚀 Quick Start

1. **Clone the repository**
   ```bash
   git clone https://github.com/yourusername/pysio_ai.git
   cd pysio_ai
   ```

2. **Set up your environment**
   ```bash
   # Install everything (Python, virtual environment, and dependencies)
   make install
   
   # Activate the virtual environment
   source env/bin/activate
   ```

3. **Configure your API key**
   ```bash
   # Create a .env file
   echo "OPENAI_API_KEY=your-api-key-here" > .env
   ```

4. **Run the application**
   ```bash
   python main.py
   ```

5. **Access the web interface**
   Open your browser and navigate to: `http://localhost:8000`

## 📁 Project Structure

```
pysio_ai/
├── api/
│   └── engine/
│       ├── diagnosis_assistant.py
│       └── evidence_retrieval.py
├── templates/
│   ├── static/
│   │   └── css/
│   │       └── styles.css
│   └── index.html
├── .env
├── main.py
├── Makefile
└── README.md
```

## 🛠 Development Commands

The project includes a Makefile with several useful commands:

```bash
make install     # Set up everything (pyenv, virtualenv, and dependencies)
make clean       # Remove virtual environment and cached files
make setup-pyenv # Install Python 3.11.5 using pyenv
make setup-venv  # Create and activate virtual environment
make install-deps# Install project dependencies
make help       # Show available commands
```

## 🔧 Configuration

The application requires an OpenAI API key to function. Create a `.env` file in the project root with:

```plaintext
OPENAI_API_KEY=your-api-key-here
```

## 🚀 Usage

1. Start the server using `python main.py`
2. Open your web browser to `http://localhost:8000`
3. Enter symptoms in the text area
4. Click "Analyze Symptoms" to receive an analysis

## 🔒 Security Notes

- Never commit your `.env` file or expose your API keys
- The application performs validation checks on startup
- All API requests are made server-side for security

## 🐛 Troubleshooting

Common issues and solutions:

1. **Server won't start**
   - Check if `.env` file exists with valid API key
   - Ensure virtual environment is activated
   - Verify all dependencies are installed

2. **Missing dependencies**
   ```bash
   make install-deps
   ```

3. **Python version mismatch**
   ```bash
   make clean
   make install
   ```

## 📝 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## 📬 Contact

Oscar Valls Lozano - (https://www.linkedin.com/in/oscar-valls-lozano/)

Project Link: [https://github.com/osscarvalls/pysio_ai](https://github.com/osscarvalls/pysio_ai)

## 📚 Entrez/PubMed Access

This application uses NCBI's E-utilities to retrieve medical literature from PubMed. While an API key is not required, it's good practice to:

1. **Provide an Email Address**
   Add this to your `.env` file:
   ```
   ENTREZ_EMAIL=your-email@example.com
   ```
   This allows NCBI to contact you if there are problems with your requests.

2. **Usage Guidelines**
   - Default limit: 3 requests/second
   - Please respect NCBI's [Usage Guidelines](https://www.ncbi.nlm.nih.gov/books/NBK25497/)
   - If you need to make more intensive requests, consider [obtaining an API key](https://www.ncbi.nlm.nih.gov/account/settings/)

Note: The application will work without these configurations, but providing an email address is considered courteous to NCBI's services.
