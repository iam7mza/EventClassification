import pandas as pd
import numpy as np
import os
import glob
import ast

def create_event_type_mapping(csv_files):
    """
    Create a dictionary mapping event types to target numbers
    
    Args:
        csv_files: List of CSV file paths
    
    Returns:
        dict: Dictionary mapping event_type to target number
    """
    event_to_target = {}
    current_target = 0
    
    # Extract unique event types from filenames
    for csv_file in csv_files:
        filename = os.path.basename(csv_file)  # Get just the filename
        event_name = filename.replace('.csv', '')  # Remove .csv extension
        # Remove _sample suffix if it exists
        if '_sample' in event_name:
            event_name = event_name.split('_sample')[0]
        event_type = event_name.split('_')[0]  # Extract the first part
        
        # Assign target if not seen before
        if event_type not in event_to_target:
            event_to_target[event_type] = current_target
            current_target += 1
    
    return event_to_target


def create_unique_event_mapping(csv_files):
    """
    Create a dictionary mapping unique particle types to target numbers
    Groups similar particles together (e.g., all Z bosons get same target)
    
    Args:
        csv_files: List of CSV file paths
    
    Returns:
        dict: Dictionary mapping event_type to target number based on unique particles
    """
    event_to_target = {}
    current_target = 0
    
    # Define particle groupings - add more as needed
    #TODO check list and update properly; also discuss with Andrea and maybe Prof#
    particle_groups = {
        'Z': ['Z', 'Zee', 'Zmumu', 'Ztautau', 'ZqqZll',],
        'W': ['W', 'Wenu', 'Wmunu', 'Wtaunu', 'Wminusenu', 'Wminusmunu', 'Wplusenu', 'Wplusmunu', 'Wplustaunu', 'WlvZqq', 'WqqZll', 'WplvWmqq', 'Wminustaunu'],
        'H': ['H', 'Higgs', 'Hbb', 'Hgg', 'Hww', 'Hzz', 'Htt'],
        'ttbar': ['ttbar', 'tt', 'ttH', 'ttW', 'ttZ'],
        'QCD': ['QCD', 'dijet', 'multijet'],
        'singletop': ['single', 'singletop', 'st', 'tchan', 'schan', 'Wt', 'wtchan'],
        'VV': ['WW', 'WZ', 'ZZ', 'VV', 'diboson', 'llvv', 'lvvv', 'lllv'],
        'TT': ['TT'],  # Direct TT production
        'GG': ['GG'],  # Gluon-gluon processes
        'C1N2': ['C1N2'],  # SUSY processes
        'ZPrime': ['ZPrime400', 'ZPrime500', 'ZPrime750', 'ZPrime1000', 'ZPrime1250', 'ZPrime1500', 'ZPrime1750', 'ZPrime2000', 'ZPrime2250','ZPrime2500', 'ZPrime2750','ZPrime3000']
    }
    
    # Create reverse mapping: specific event -> particle group
    event_to_particle = {}
    for particle, variations in particle_groups.items():
        for variation in variations:
            event_to_particle[variation] = particle
    
    # Extract unique particle types from filenames
    unique_particles = set()
    for csv_file in csv_files:
        filename = os.path.basename(csv_file)
        event_name = filename.replace('.csv', '')  # Remove .csv extension
        # Remove _sample suffix if it exists
        if '_sample' in event_name:
            event_name = event_name.split('_sample')[0]
        
        event_type = event_name.split('_')[0]
        
        # Find which particle group this belongs to
        particle = event_to_particle.get(event_type, event_type)  # Default to itself if not found
        unique_particles.add(particle)
    
    # Assign targets to unique particles
    particle_to_target = {}
    for particle in sorted(unique_particles):
        particle_to_target[particle] = current_target
        current_target += 1
    
    # Create final mapping for all event types
    for csv_file in csv_files:
        filename = os.path.basename(csv_file)
        event_name = filename.replace('.csv', '')  # Remove .csv extension
        # Remove _sample suffix if it exists
        if '_sample' in event_name:
            event_name = event_name.split('_sample')[0]
            
        event_type = event_name.split('_')[0]
        
        # Map to particle group, then to target
        particle = event_to_particle.get(event_type, event_type)
        event_to_target[event_type] = particle_to_target[particle]
    
    return event_to_target


def target_to_vector(targetDict):
    """
    Convert a dictionary mapping event types to target numbers into one-hot encoded vectors
    
    Args:
        targetDict: Dictionary mapping event_type to target number
    
    Returns:
        dict: Dictionary mapping event_type to one-hot encoded vector
    """
    num_classes = len(set(targetDict.values()))  # Number of unique targets
    target_to_onehot = {}
    
    for event_type, target in targetDict.items():
        one_hot = np.zeros(num_classes)
        one_hot[target] = 1
        target_to_onehot[event_type] = one_hot
    
    return target_to_onehot



# csv_files = glob.glob('csv_output/*.csv')

# targets = create_unique_event_mapping(csv_files)
# targets = target_to_vector(targets) #converting to one-hot encoded vectors


##TODO : create a function that make the df with the targets###
####### should be similar to the event_classification.py file ########

def load(csv_files_path, unique=True, targettype='onehot'):
    """
    Loads CSV files and returns a dataframe with all event and corresponding targets
    Args:
        csv_files_path: Path to CSV files
        unique: If True, uses unique event mapping (i.e. all events with initial state Z will have the same target); otherwise uses event type mapping
        targettype: Type of target to return; 'onehot' for one-hot encoded vectors, 'num' for numeric targets

    Returns:
        pd.DataFrame: DataFrame containing all events and their targets, 
    """
    
    csv_files = glob.glob(f'{csv_files_path}/*.csv')
    if not csv_files:
        print("No CSV files found in csv_output directory")
        return None
    else:
        print(f"Found {len(csv_files)} CSV files")

    if unique:
        targets_numeric = create_unique_event_mapping(csv_files)
    else:
        targets_numeric = create_event_type_mapping(csv_files)

    # Keep numeric targets for assignment to DataFrame
    if targettype == 'onehot':
        # We'll create one-hot vectors later, but use numeric for now
        targets_onehot = target_to_vector(targets_numeric)
        print(f"One-hot encoded targets created with {len(set(targets_numeric.values()))} classes")
    elif targettype == 'num':
        targets_onehot = None
    else:
        raise ValueError("Invalid targettype. Use 'onehot' or 'num'.")

    print(f"Event type mapping:")
    for event_type, target in sorted(targets_numeric.items()):
        print(f"  {event_type}: target {target}")

    # List to store all dataframes
    dataframes = []

    # Process each CSV file using the mapping
    for csv_file in csv_files:
        print(f"Processing: {csv_file}")
        
        # Read the CSV file
        df = pd.read_csv(csv_file)
        
        # Define all list columns based on the new CSV structure
        list_columns = [
            # Lepton vectors
            'lep_truthMatched', 'lep_trigMatched', 'lep_pt', 'lep_eta', 'lep_phi', 'lep_E', 
            'lep_z0', 'lep_charge', 'lep_type', 'lep_isTightID', 'lep_ptcone30', 'lep_etcone20', 
            'lep_trackd0pvunbiased', 'lep_tracksigd0pvunbiased',
            
            # Jet vectors
            'jet_pt', 'jet_eta', 'jet_phi', 'jet_E', 'jet_jvt', 'jet_trueflav', 
            'jet_truthMatched', 'jet_MV2c10',
            
            # Photon vectors
            'photon_truthMatched', 'photon_trigMatched', 'photon_pt', 'photon_eta', 'photon_phi', 
            'photon_E', 'photon_isTightID', 'photon_ptcone30', 'photon_etcone20', 'photon_convType',
            
            # Large R jet vectors
            'largeRjet_pt', 'largeRjet_eta', 'largeRjet_phi', 'largeRjet_E', 'largeRjet_m', 
            'largeRjet_truthMatched', 'largeRjet_D2', 'largeRjet_tau32',
            
            # Tau vectors
            'tau_pt', 'tau_eta', 'tau_phi', 'tau_E', 'tau_charge', 'tau_isTightID', 
            'tau_truthMatched', 'tau_trigMatched', 'tau_nTracks', 'tau_BDTid',
            
            # Systematic uncertainty vectors
            'lep_pt_syst', 'jet_pt_syst', 'photon_pt_syst', 'largeRjet_pt_syst', 'tau_pt_syst'
        ]
        
        # Convert string representations of lists to actual lists
        for col in list_columns:
            if col in df.columns:
                df[col] = df[col].apply(ast.literal_eval)
            else:
                print(f"Warning: Column {col} not found in CSV")
        
        # Extract event name from filename
        filename = os.path.basename(csv_file)
        event_name = filename.replace('.csv', '')  # Remove .csv extension
        # Remove _sample suffix if it exists
        if '_sample' in event_name:
            event_name = event_name.split('_sample')[0]
            
        event_type = event_name.split('_')[0]
        
        # Get numeric target for DataFrame assignment
        target_numeric = targets_numeric[event_type]
        
        # Add columns to dataframe
        df['target'] = target_numeric  # Use numeric target for all rows
        df['event_type'] = event_type
        df['full_event_name'] = event_name
        
        # If one-hot encoding requested, add one-hot vector as a separate column
        if targettype == 'onehot':
            onehot_vector = targets_onehot[event_type]
            df['target_onehot'] = [onehot_vector] * len(df)  # Same vector for all rows in this file
        
        # Append to list
        dataframes.append(df)
        print(f"  -> Added {len(df)} events with target {target_numeric} (event type: {event_type})")

    # Combine all dataframes
    if dataframes:
        combined_df = pd.concat(dataframes, ignore_index=True)
        print(f"\nCombined dataset: {len(combined_df)} total events")
        
        print(f"\nTarget distribution:")
        print(combined_df['target'].value_counts().sort_index())
        
        # Display first few rows with updated column names
        print(f"\nFirst few rows:")
        if targettype == 'onehot':
            display_cols = ['event', 'runNumber', 'eventNumber', 'channelNumber', 'jet_n', 'lep_n', 
                          'target', 'target_onehot', 'event_type', 'full_event_name']
        else:
            display_cols = ['event', 'runNumber', 'eventNumber', 'channelNumber', 'jet_n', 'lep_n', 
                          'target', 'event_type', 'full_event_name']
        
        # Only show columns that exist in the dataframe
        available_cols = [col for col in display_cols if col in combined_df.columns]
        print(combined_df[available_cols].head())
        
    else:
        print("No data loaded")
        return None


    print(f"\nFinal dataset shape: {combined_df.shape}")
    return combined_df


def clean(df):
    """
    Cleans dataframe by removing any values that are empty list or lists with only one value
    [] -> [0,0] and [val] -> [val,0]

    args:
        df: DataFrame to clean
    returns:
        pd.DataFrame: Cleaned DataFrame
    """
    list_columns = []
    for col in df.columns:
        if df[col].dtype == 'object':  # String columns that might contain lists
            sample_values = df[col].dropna().head(5).values
            if any('[' in str(val) for val in sample_values):
                list_columns.append(col)

    if list_columns == []:
        print("No list columns found to clean.")
        return df

    def fix_list_cell(x):
        # parse the string into a Python list if needed
        if isinstance(x, str):
            try:
                x = ast.literal_eval(x)
            except (ValueError, SyntaxError):
                return x
        # now x is either a list or something else
        if isinstance(x, list):
            if len(x) == 0:
                return [0, 0]
            if len(x) == 1:
                return [x[0], 0]
        return x

    for col in list_columns:
        df[col] = df[col].apply(fix_list_cell)

    print(f"Cleaned {len(list_columns)} list columns in DataFrame")
    return df