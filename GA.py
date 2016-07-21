def Preprocess():
    '''text'''
    
def Prediction():
    '''text'''
    
def GeneticAlgorithm():
    '''text'''

def main(args):
    #Process the Training Data
    Preprocess(args.Training)
    #Define the algorithm
    
    #Predict with the Scoring Data


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Search for the best prediction model using a Genetic algorithm')
    parser.add_argument('Training',help='The file that the training will be based on')
    parser.add_argument('Scoring',help='The values that are needed to predict')
    parser.add_argument('Prediction_method',choices=[''],help='The prediction method to use')
    args = parser.parse_args()
    print(args)
