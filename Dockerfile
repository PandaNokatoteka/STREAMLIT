# Βασική εικόνα με Python
FROM python:3.10-slim

# Ορισμός μεταβλητών περιβάλλοντος
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# Εγκατάσταση απαραίτητων εξαρτήσεων συστήματος
RUN apt-get update && apt-get install -y     build-essential     python3-dev     libhdf5-dev     libz-dev     libglib2.0-0     libsm6     libxext6     libxrender-dev     && rm -rf /var/lib/apt/lists/*

# Δημιουργία φακέλου εργασίας
WORKDIR /app

# Αντιγραφή αρχείων εφαρμογής
COPY . /app

# Εγκατάσταση Python εξαρτήσεων
RUN pip install --upgrade pip && pip install -r requirements.txt

# Έναρξη της εφαρμογής
CMD ["streamlit", "run", "streamlit_app.py", "--server.port=8501", "--server.enableCORS=false"]