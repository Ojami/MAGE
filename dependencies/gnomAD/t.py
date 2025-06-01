#!/usr/bin/env python

import requests
import pprint
prettyprint = pprint.PrettyPrinter(indent=2).pprint

def fetch(jsondata, url="https://gnomad.broadinstitute.org/api"):
    # The server gives a generic error message if the content type isn't
    # explicitly set
    headers = {"Content-Type": "application/json"}
    response = requests.post(url, json=jsondata, headers=headers)
    json = response.json()
    if "errors" in json:
        raise Exception(str(json["errors"]))
    return json

def get_variant_list(gene_id, dataset="gnomad_r2_1"):
    # Note that this is GraphQL, not JSON.
    fmt_graphql = """
    {
      variant(variantId: "%s", dataset: %s) {
        variantId
        sortedTranscriptConsequences {
          lof
        }
      }
    }
    """
    # This part will be JSON encoded, but with the GraphQL part left as a
    # glob of text.
    req_variantlist = {
        "query": fmt_graphql % (gene_id, dataset),
        "variables": {}
        }
    #response = fetch(req_variantlist)
    return req_variantlist

prettyprint(get_variant_list("1-55517991-C-CAT"))


