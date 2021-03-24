function varargout = Parity_Checker(varargin)
% PARITY_CHECKER MATLAB code for Parity_Checker.fig
%      PARITY_CHECKER, by itself, creates a new PARITY_CHECKER or raises the existing
%      singleton*.
%
%      H = PARITY_CHECKER returns the handle to a new PARITY_CHECKER or the handle to
%      the existing singleton*.
%
%      PARITY_CHECKER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PARITY_CHECKER.M with the given input arguments.
%
%      PARITY_CHECKER('Property','Value',...) creates a new PARITY_CHECKER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Parity_Checker_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Parity_Checker_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Parity_Checker

% Last Modified by GUIDE v2.5 25-May-2013 17:22:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Parity_Checker_OpeningFcn, ...
                   'gui_OutputFcn',  @Parity_Checker_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Parity_Checker is made visible.
function Parity_Checker_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Parity_Checker (see VARARGIN)

% Choose default command line output for Parity_Checker
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Parity_Checker wait for user response (see UIRESUME)
 uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Parity_Checker_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
B = get(hObject, 'Data');

assignin('base','B',B)



% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Parity_Generator(B);




































function [ C ] = Parity_Generator( A )
%Parity_Generator Summary of this function goes here
%   Detailed explanation goes here

for I =1:4
    
    for J = 1:4
        temp = temp + A(I,J);
    end
    if (mod(temp,2) ~= 1)
            A(I,5)=0;
    else
            A(I,5)=1;
    end
    temp = 0;
end
for I =1:4
    temp = 0;
    for J = 1:4
        temp = temp + A(J,I);
    end
    if (mod(temp,2)~= 1)
            A(5,I)=0;
    else
            A(5,I)=1;
    end
end
C = A;

