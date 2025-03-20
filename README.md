# Fitness Assistant Application

## Overview
This application is a **Fitness Assistant** that combines real-time workout tracking using computer vision with scientific research integration to provide feedback and generate a report on the effectiveness of specific exercises (e.g., bicep curls) in muscle hypertrophy. It uses **MediaPipe** for pose estimation, **PubMed** for fetching relevant scientific articles, and **Hugging Face Transformers** for generating a report based on workout data and research insights.

---

## Features
1. **Real-time Workout Tracking**:
   - Tracks bicep curl exercises using a webcam.
   - Counts repetitions (reps) and sets.
   - Measures **time under tension** for concentric (lifting) and eccentric (lowering) phases.
   - Records the **range of motion** (min and max angles) for each rep.

2. **Scientific Research Integration**:
   - Fetches relevant articles from **PubMed** using the **Entrez API**.
   - Searches for articles related to:
     - Bicep curls.
     - Exercise angles.
     - Time under tension.
     - Muscle hypertrophy.

3. **Automated Report Generation**:
   - Combines workout data with fetched PubMed articles.
   - Uses a **Hugging Face text-generation model** (e.g., GPT-2) to generate a detailed report on the effectiveness of the workout.

4. **Real-time Feedback**:
   - Displays rep count, set count, and feedback on the screen.
   - Visualizes body landmarks and angles using **MediaPipe**.

---

## Requirements
To run this application, you need the following Python libraries installed:
- `opencv-python`
- `mediapipe`
- `numpy`
- `biopython`
- `transformers`
- `torch` (for Hugging Face models)

You can install the required libraries using:
```bash
pip install opencv-python mediapipe numpy biopython transformers torch
