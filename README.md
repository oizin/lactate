# Background

This repo contains code for the paper: The relationship between hyperglycaemia on admission and patient outcome is modified by hyperlactatemia and diabetic status: a retrospective analysis of the eICU collaborative research database.

# Using this repository

## eICU

As we are not data custodians we cannot publicaly share the eICU data used. However, it is available to those with credentialed access to physionet.org. Credentialed access can be requested through your physionet account.

We used dbt to connect to the Google Bigquery eICU database that can be autopopulated through the physionet interface.

https://docs.getdbt.com/reference/warehouse-setups/bigquery-setup

## Data extraction  and nalysis

The analysis uses R, it can be repeated by:

```
./run.sh
```