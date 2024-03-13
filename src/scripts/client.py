import urllib.request
import json
import requests
import csv
import time

url_tissue_rui = "https://cdn.humanatlas.io/digital-objects/ds-graph/hra-pop-full/v0.5.1/assets/full-dataset-graph.jsonld"
url_reference_organ_data = "https://hubmapconsortium.github.io/ccf-ui/rui/assets/reference-organ-data.json"
url_collision_detection = "http://0.0.0.0:8080/get-collisions"


t1 = time.time()
fileds = ['rui_location_id', 'number_of_collisions', 'CPU_time', 'GPU_time']

rows = []
with urllib.request.urlopen(url_tissue_rui) as content:
    data = json.loads(content.read().decode())
    graph = data['@graph']

    
    for person in graph:
        samples = person['samples']
        for sample in samples:
            if sample['sample_type'] == 'Tissue Block':
                if 'rui_location' not in sample:
                    print(sample['@id'],"doesn't have rui location-----------------------")
                    continue

                rui_location = sample['rui_location']
                rui_location_id = rui_location['@id']
                
                r = requests.post(url=url_collision_detection, json=rui_location)
                if r:
                    response = r.json()
                    print(rui_location_id)
                    if 'number_of_collisions' not in response:
                        print(response)
                        continue
                    number_of_collisions = response['number_of_collisions']
                    print('The rui location collides with ', number_of_collisions)

                    if 'CPU_time' in response:
                        CPU_time = response['CPU_time']
                        GPU_time = response['GPU_time']
                        print('CPU_time: ', CPU_time, ', GPU_time: ', GPU_time)
                        row = [rui_location_id, number_of_collisions, CPU_time, GPU_time]
                    else:
                        row = [rui_location_id, number_of_collisions, '-', '-']
                        
                    rows.append(row[:])
        
        break
        


with open('report.txt', 'w') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(fileds)
    writer.writerows(rows) 

t2 = time.time()

print("total running time:", t2 - t1)
