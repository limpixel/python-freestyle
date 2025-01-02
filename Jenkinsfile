pipeline {
    agent any
    stages {
        stage('Build Docker Image') {
            steps {
                script {
                    def branchName = env.BRANCH_NAME
                    echo "Building Docker image for branch: ${branchName}"
                    bat "docker build -t python-app:%BRANCH_NAME% ."
                }
            }
        }
        stage('Run Docker Container') {
            steps {
                script {
                    echo "Running Docker container..."
                    bat "docker run --rm python-app:%BRANCH_NAME%"
                }
            }
        }
    }
}