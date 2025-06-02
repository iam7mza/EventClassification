import pandas as pd
import numpy as np
import os
import glob
import ast

def load(csv_files_path):
    """
    Loads CSV files and returns a dataframe with all event and corresponding targets
    Args:
        csv_files_path: Path to CSV files

    Returns:
        pd.DataFrame: DataFrame containing all events and their targets, 
    """
    
    csv_files = glob.glob(f'{csv_files_path}/*.csv')
    if not csv_files:
        print("No CSV files found in csv_output directory")
        return None
    else:
        print(f"Found {len(csv_files)} CSV files")

    # if unique:
    #     targets_numeric = create_unique_event_mapping(csv_files)
    # else:
    #     targets_numeric = create_event_type_mapping(csv_files)

    # Keep numeric targets for assignment to DataFrame
    # if targettype == 'onehot':
    #     # We'll create one-hot vectors later, but use numeric for now
    #     targets_onehot = target_to_vector(targets_numeric)
    #     print(f"One-hot encoded targets created with {len(set(targets_numeric.values()))} classes")
    # elif targettype == 'num':
    #     targets_onehot = None
    # else:
    #     raise ValueError("Invalid targettype. Use 'onehot' or 'num'.")

    # print(f"Event type mapping:")
    # for event_type, target in sorted(targets_numeric.items()):
    #     print(f"  {event_type}: target {target}")

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
            
    
        # Add column to dataframe
        df['full_event_name'] = event_name
        
        # Append to list
        dataframes.append(df)
        print(f"  -> Added {len(df)} events (event name: {event_name})")

    # Combine all dataframes
    if dataframes:
        combined_df = pd.concat(dataframes, ignore_index=True)
        print(f"\nCombined dataset: {len(combined_df)} total events")
        print(f"Columns: {combined_df.columns.tolist()}")
        # Display first few rows with updated column names
        print(f"\nFirst few rows:")
        display_cols = ['event', 'runNumber', 'eventNumber', 'channelNumber', 'jet_n', 'lep_n', 'full_event_name']

        
        # Only show columns that exist in the dataframe
        available_cols = [col for col in display_cols if col in combined_df.columns]
        print(combined_df[available_cols].head())
        
    else:
        print("No data loaded")
        return None


    print(f"\nFinal dataset shape: {combined_df.shape}")

    # Extract target values
    combined_df = target_extract(combined_df)
    return combined_df



def target_extract(df):
    """
    Extracts target values from the DataFrame based on the "channelNumber" column.
    Args:
        df (pd.DataFrame): DataFrame containing the data.
    Returns:
        df (pd.DataFrame): DataFrame with additional target columns.
    """

    # Define particle grouips (channel numbers corresonding to which event)
    #got channel number from https://www.pi.uni-bonn.de/brock/en/results/data/t00000029.pdf p65
    particle_groups = {}

    # ttbar production
    particle_groups[410000] = 'ttbar'

    # single top production (all variants)
    single_top_channels = [410011, 410012, 410013, 410014, 410025, 410026, 410062, 410063, 410050]
    for channel in single_top_channels:
        particle_groups[channel] = 'singletop'

    # Z boson production (original range)
    for channel in [361100, 361101, 361102, 361103, 361104, 361105, 361106, 361107, 361108]:
        particle_groups[channel] = 'Z'

    # Z+jets production (original range)
    for i in range(361400, 361442):
        particle_groups[i] = 'Z'

    # Z boson in association with jets (from table: 364100-364141)
    for i in range(364100, 364142):
        particle_groups[i] = 'Z'

    # W boson production (original range)
    for i in range(361500, 361506):
        particle_groups[i] = 'W'

    # W boson in association with jets (from table: 364156-364197)
    for i in range(364156, 364198):
        particle_groups[i] = 'W'

    # Diboson production (original)
    for channel in [363359, 363360, 363492, 363356, 363490, 363358, 363489, 363491, 363493]:
        particle_groups[channel] = 'VV'

    # Diboson production (from table)
    diboson_from_table = [361600, 361601, 361602, 361603, 361604, 361606, 361607, 361609, 361610]
    for channel in diboson_from_table:
        particle_groups[channel] = 'VV'

    # Higgs production
    for channel in [345324, 345323, 345060, 344235, 341947, 341964, 343981, 345041, 345318, 345319, 341081]:
        particle_groups[channel] = 'H'

    # BSM production
    particle_groups[301325] = 'ZPrime'
    particle_groups[392985] = 'SUSY'

    #C1N2
    particle_groups[392220] = 'C1N2'
    particle_groups[392217] = 'C1N2'
    particle_groups[392223] = 'C1N2'
    
    #gluon gluon
    particle_groups[370114] = 'GG'
    particle_groups[370118] = 'GG'
    particle_groups[370129] = 'GG'
    particle_groups[370144] = 'GG'

    #Zprime750
    particle_groups[301324] = 'Zprime750'

    #Zprime2000
    particle_groups[301329] = 'Zprime2000'


    # Add target columns based on channelNumber
    channels = np.unique(df['channelNumber'])
    uniqueTargets = np.unique([particle_groups[i] for i in channels if i in particle_groups])
    nTargets = len(uniqueTargets)

    # Apply the function to get all three return values
    def process_channel(channel):
        return channel_to_onehot(channel, particle_groups, uniqueTargets)
    
    # Get results for all channels
    results = df['channelNumber'].apply(process_channel)
    
    # Extract the three components
    df['onehot_target'] = results.apply(lambda x: x[0])
    df['numeric_target'] = results.apply(lambda x: x[1])
    df['event_type'] = results.apply(lambda x: x[2])
    
    print(f"\nTarget extraction complete:")
    print(f"Unique particle types: {uniqueTargets}")
    print(f"Number of classes: {nTargets}")
    
    # Show distribution
    print(f"\nParticle type distribution:")
    print(df['event_type'].value_counts())
    
    return df




def channel_to_onehot(channel_number, particle_groups, uniqueTargets):
    """
    Convert channel number to one-hot encoded target vector, numerical target, and particle type
    
    Args:
        channel_number: The channel number (or Series) to convert
        particle_groups: Dictionary mapping channel numbers to particle types
        uniqueTargets: Array of unique targets found in the actual data
    
    Returns:
        tuple: (onehot_vector, numerical_target, particle_type) or (None, -1, 'Unknown') if channel not found
    """
    import numpy as np
    import pandas as pd
    
    # Handle Series input (for apply operations)
    if isinstance(channel_number, pd.Series):
        return channel_number.apply(lambda x: channel_to_onehot(x, particle_groups, uniqueTargets))
    
    # Use the uniqueTargets from your data, not all possible particle types
    particle_to_index = {particle: idx for idx, particle in enumerate(sorted(uniqueTargets))}
    
    # Check if channel exists in particle_groups
    if channel_number not in particle_groups:
        print(f"Warning: Channel {channel_number} not found in particle_groups")
        return None, -1, 'Unknown'
    
    # Get particle type for this channel
    particle_type = particle_groups[channel_number]
    
    # Check if this particle type is in our actual data
    if particle_type not in particle_to_index:
        print(f"Warning: Particle type {particle_type} not found in uniqueTargets")
        return None, -1, particle_type
    
    # Get index for this particle type
    target_index = particle_to_index[particle_type]
    
    # Create one-hot vector with correct size
    onehot = np.zeros(len(uniqueTargets))
    onehot[target_index] = 1
    
    return onehot, target_index, particle_type





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
