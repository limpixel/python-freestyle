pipeline {
    agent any
    environment {
        BRANCH_NAME = "${env.BRANCH_NAME ?: 'latest'}"
    }
    stages {

        stage('checking branch') {
            steps {
                script {
                    echo "Current branch: ${BRANCH_NAME}"
                }
            }
        }

        stage('Build Docker Image') {
            steps {
                script {
                    echo "Building Docker image for branch: ${BRANCH_NAME}"
                    sh "docker build -t python-app:${BRANCH_NAME} ."
                }
            }
        }
        stage('Run Docker Container') {
            steps {
                script {
                    echo "Running Docker container..."
                    sh "docker run --rm python-app:${BRANCH_NAME}"
                }
            }
        }
    }
}
