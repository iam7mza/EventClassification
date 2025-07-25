# EventClassification
Event Classification using ML on Atlas Open Data sample (Still under development)

## Dataset
The dataset used for this project is available at: [Atlas Open Data](https://atlas-opendata.web.cern.ch/Legacy13TeV/1largeRjet1lep/).

## Files

- `EventClassification.ipynb`: Contains the code for training the models.
- `CorrelationTest.ipynb`: Used for testing correlations in the data.
- `rootToCSV`: No one knows what this one does.

The *1 lepton + 1 large jet* data is located in the `/Data` directory. Additional *2 lepton* data is available in the `2lepData` directory.

## How to Use `rootToCSV`

To convert a `.root` file to `.csv`, run the following command:
```bash
./rootToCSV [file_name.root] [number_of_samples]
```

To process all `.root` files in a directory, execute:
```bash
source script.sh
```
This command should be run in the directory containing the `.root` files. You can specify the maximum number of samples to extract from each file by modifying the `SAMPLE_SIZE` variable in `script.sh`. The resulting `.csv` files will be moved to the `/csv_output` directory.

**NB:** due to an unresoved issue not all `.csv` files may be moved to the `/csv_output` directory. In such case run
```bash
mv *.csv csv_output/.
```

# License:

This model is not good *(yet)*. We suggest you find a different one to use.