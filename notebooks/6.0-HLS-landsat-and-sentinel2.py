import requests 

response = requests.post('https://appeears.earthdatacloud.nasa.gov/api/login', auth=('', '!'))
token_response = response.json()
print(token_response)

token = token_response['token']
task_id_list = ['16760cff-b2a0-497e-bbed-64d1f35387eb', #landsat Part 1
                '7d9a8aca-27c9-44f0-9a8f-4b4921c2fd64', #landsat Part 2
                '19cfbafc-335e-40f4-adf8-10ae72ba42fb', #landsat Part 3
                'fb0daa9e-c6bc-4223-afd3-030c1627ab70', #landsat Part 4
                '04967a71-498b-412e-82e1-929ae0d623ea', #Sentinel 2a/2b Part 1
                'c9182642-49be-41bf-8100-abc5381480b8', #Sentinel 2a/2b Part 2
                'ea0339a3-966b-49f9-8289-ec8891a6d254', #Sentinel 2a/2b Part 3
                '76e59a18-ca26-47d5-b89c-3f06aa9a8939', #Sentinel 2a/2b Part 4
                'bb3297dc-be0c-4b07-9b7a-ef646f5f1b89', #Sentinel 2a/2b Part 5
                '39c9be31-5715-4b96-9774-847c27c731f2', #Sentinel 2a/2b Part 6
                '4d09fbad-fc53-4644-829c-48fa023147e0', #Sentinel 2a/2b Part 7
                '5cbb5b7a-b714-4ef4-8170-a61302e1a6ea'] #Sentinel 2a/2b Part 8

import os
token = token_response['token']

for j in range(len(task_id_list)):
    task_id = task_id_list[j]
    response = requests.get(
        'https://appeears.earthdatacloud.nasa.gov/api/bundle/{0}'.format(task_id),  
        headers={'Authorization': 'Bearer {0}'.format(token)}
    )
    bundle_response = response.json()
    for i in range(len(bundle_response['files'])):
        print(i)
        # get a stream to the bundle file
        file_id = bundle_response['files'][i]['file_id']
        filename = bundle_response['files'][i]['file_name']
        response = requests.get( 
            'https://appeears.earthdatacloud.nasa.gov/api/bundle/{0}/{1}'.format(task_id,file_id),  
            headers={'Authorization': 'Bearer {0}'.format(token)}, 
            allow_redirects=True,
            stream=True
        ) 

        # create a destination directory to store the file in
        # dest_dir = "C:\\Users\\admin\\Desktop\\GitHub\\riskaudit\\data\\raw\\HLS_landsat" #Landsat
        # dest_dir = "C:\\Users\\admin\\Desktop\\GitHub\\riskaudit\\data\\raw\\HLS_sentinel2" #Sentinel 2a/2b
        filepath = os.path.join(dest_dir, filename)
        os.makedirs(os.path.dirname(filepath), exist_ok=True)

        # write the file to the destination directory
        with open(filepath, 'wb') as f:
            for data in response.iter_content(chunk_size=8192):
                f.write(data)