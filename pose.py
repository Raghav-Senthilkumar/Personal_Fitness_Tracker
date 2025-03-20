import json
import time
import cv2
import numpy as np
import mediapipe as mp
from Bio import Entrez
from transformers import pipeline

Entrez.email = "your_email@example.com"

query = "bicep curls OR exercise angles OR time under tension OR muscle hypertrophy"
handle = Entrez.esearch(db="pubmed", term=query, retmax=100)
record = Entrez.read(handle)
handle.close()

pmids = record["IdList"]
print("PubMed IDs:", pmids)

def fetch_article_details(pmids):
    handle = Entrez.efetch(db="pubmed", id=pmids, rettype="abstract", retmode="text")
    articles = handle.read()
    handle.close()
    return articles

mp_pose = mp.solutions.pose
pose = mp_pose.Pose(static_image_mode=False, min_detection_confidence=0.5, min_tracking_confidence=0.5)
mp_drawing = mp.solutions.drawing_utils

def calculate_angle(a, b, c):
    a = np.array(a)
    b = np.array(b)
    c = np.array(c)
    radians = np.arctan2(c[1] - b[1], c[0] - b[0]) - np.arctan2(a[1] - b[1], a[0] - b[0])
    angle = np.abs(radians * 180.0 / np.pi)
    if angle > 180.0:
        angle = 360 - angle
    return angle

cap = cv2.VideoCapture(0)

exercise_type = "bicep_curl"
rep_count = 0
set_count = 0
stage = None
feedback = "Start your workout!"
workout_data = {
    "exercise_type": exercise_type,
    "sets": []
}
current_set = {
    "reps": []
}
start_time = None
min_angle = 180
max_angle = 0
concentric_time = 0
eccentric_time = 0

while cap.isOpened():
    ret, frame = cap.read()
    if not ret:
        break

    image = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
    image.flags.writeable = False
    results = pose.process(image)
    image.flags.writeable = True
    image = cv2.cvtColor(image, cv2.COLOR_RGB2BGR)

    try:
        landmarks = results.pose_landmarks.landmark
        shoulder = [landmarks[mp_pose.PoseLandmark.LEFT_SHOULDER.value].x,
                    landmarks[mp_pose.PoseLandmark.LEFT_SHOULDER.value].y]
        elbow = [landmarks[mp_pose.PoseLandmark.LEFT_ELBOW.value].x,
                 landmarks[mp_pose.PoseLandmark.LEFT_ELBOW.value].y]
        wrist = [landmarks[mp_pose.PoseLandmark.LEFT_WRIST.value].x,
                 landmarks[mp_pose.PoseLandmark.LEFT_WRIST.value].y]
        angle = calculate_angle(shoulder, elbow, wrist)
        min_angle = min(min_angle, angle)
        max_angle = max(max_angle, angle)

        if angle > 120 and stage != "down":
            if stage == "up":
                concentric_time = time.time() - start_time if start_time else 0
            stage = "down"
            start_time = time.time()
        elif angle < 60 and stage == "down":
            stage = "up"
            eccentric_time = time.time() - start_time if start_time else 0
            rep_count += 1
            current_set["reps"].append({
                "rep_number": rep_count,
                "time_concentric": concentric_time,
                "time_eccentric": eccentric_time,
                "angle_of_full_flex": min_angle,
                "angle_of_full_retraction": max_angle
            })
            start_time = time.time()
            min_angle = 180
            max_angle = 0

            if rep_count == 8:
                set_count += 1
                workout_data["sets"].append(current_set)
                rep_count = 0
                current_set = {
                    "reps": []
                }

                if set_count == 3:
                    print(json.dumps(workout_data, indent=4))
                    break

    except:
        pass

    cv2.putText(image, feedback, (10, 30), cv2.FONT_HERSHEY_SIMPLEX, 1, (0, 255, 0), 2, cv2.LINE_AA)
    cv2.putText(image, f"Reps: {rep_count}", (10, 70), cv2.FONT_HERSHEY_SIMPLEX, 1, (0, 255, 0), 2, cv2.LINE_AA)
    cv2.putText(image, f"Set: {set_count + 1}/3", (10, 110), cv2.FONT_HERSHEY_SIMPLEX, 1, (0, 255, 0), 2, cv2.LINE_AA)
    mp_drawing.draw_landmarks(image, results.pose_landmarks, mp_pose.POSE_CONNECTIONS)
    cv2.imshow("Fitness Assistant", image)

    if cv2.waitKey(10) & 0xFF == ord('q'):
        break

cap.release()
cv2.destroyAllWindows()

articles = fetch_article_details(pmids)

report_generator = pipeline('text-generation', model='gpt2')

prompt = f"Workout Data:\n{json.dumps(workout_data, indent=4)}\n\nArticles:\n{articles}\n\nGenerate a report on the effectiveness of bicep curls, exercise angles, and time under tension in muscle hypertrophy."

report = report_generator(prompt, max_length=500, num_return_sequences=1)[0]['generated_text']
print("\nGenerated Report:\n", report)