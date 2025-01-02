pipeline {
    agent any
    stages {
        stage('Build Docker Image') {
            steps {
                script {
                    def branchName = env.BRANCH_NAME
                    echo "Building Docker image for branch: ${branchName}"
                    sh "docker build -t python-app:%BRANCH_NAME% ."
                }
            }
        }
        stage('Run Docker Container') {
            steps {
                script {
                    echo "Running Docker container..."
                    sh "docker run --rm python-app:%BRANCH_NAME%"
                }
            }
        }
    }
}