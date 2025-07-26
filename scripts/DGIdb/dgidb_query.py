import requests
import json
import time

# Define your GraphQL query
query = """
{
  genes(names: [GENES]) {
    nodes {
      name
      interactions {
        drug {
          name
        }
        interactionScore
        interactionTypes {
          directionality
        }
        publications {
          pmid
        }
      }
    }
  }
}
"""

with open("gene_list.txt", "r") as file:
    gene_list = [line.strip() for line in file if line.strip()]
query = query.replace("GENES", ", ".join(f'"{gene}"' for gene in gene_list))

# Send the request
response = requests.post(
    'https://dgidb.org/api/graphql',
    headers={'Content-Type': 'application/json'},
    json={'query': query}
)

data = response.json()
if "errors" in data:
  print("GraphQL errors returned:")
  for error in data["errors"]:
    print(f"- {error['message']}")

timestamp = str(int(time.time()))
with open("dgidb_query_response_" + timestamp + ".js", "w") as f:
    f.write("var gene_drug_interaction = ")
    json.dump(data, f, indent=2)
