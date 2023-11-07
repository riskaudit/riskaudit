import requests

response = requests.post('https://appeears.earthdatacloud.nasa.gov/api/login', auth=('jtd33', 'Sct18627Hsc5689!'))
token_response = response.json()
print(token_response)

token = token_response['token']
task_id = '16760cff-b2a0-497e-bbed-64d1f35387eb' 
response = requests.get(
    'https://appeears.earthdatacloud.nasa.gov/api/bundle/{0}'.format(task_id),  
    headers={'Authorization': 'Bearer {0}'.format(token)}
)
bundle_response = response.json()

import os
token = token_response['token']

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
    dest_dir = "C:\\Users\\admin\\Desktop\\GitHub\\riskaudit\\data\\raw\\HLS_landsat"
    filepath = os.path.join(dest_dir, filename)
    os.makedirs(os.path.dirname(filepath), exist_ok=True)

    # write the file to the destination directory
    with open(filepath, 'wb') as f:
        for data in response.iter_content(chunk_size=8192):
            f.write(data)