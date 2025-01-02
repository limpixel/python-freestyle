# Use the official Python image from Docker Hub
FROM python:3.11-slim

# Set the working directory inside the container
WORKDIR /app

# Copy the Python application file into the container
# Salin file Python ke direktori kerja
COPY "Method-Numeric-Calculator.py" .
COPY "input.txt" .

# Install Flask (or other dependencies if necessary)
RUN pip install Flask

# Run the Python application
CMD ["python", "Method-Numeric-Calculator.py", "input.txt"]